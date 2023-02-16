/*!
 *
 * \file do_beach_estimate_actual_erosion.cpp
 * \brief Estimates actual (supply-limited) beach erosion on polygons, and constructs a between-polygon budget of actual beach sediment movement
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo
 * \date 2023
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
//#include <assert.h>

#include <cmath>
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

// /*===============================================================================================================================
//
//  Since we know the potential erosion for this polygon, this routine is able to determine actual (supply-limited) erosion on the polygon in three sediment size categories. But do not do any actual erosion yet
//
// ===============================================================================================================================*/
// int CSimulation::nEstimateBeachErosionOnPolygon(int const nCoast, int const nPoly, double const dPotentialErosionOnPolygon)
// {
//    // Totals for this polygon
//    double
//       dTotFineEroded = 0,
//       dTotSandEroded = 0,
//       dTotCoarseEroded = 0;
//
//    // Estimate how much we can erode on this polygon: traverse the polygon's shoreline in a down-coast direction (i.e. in the direction of increasing coastpoint indices), with profiles which are parallel to the polygon's up-coast boundary
//    // NOTE also tried doing a second traverse, going up-coast, with profiles parallel to the down-coast boundary, if the down-coast routine was unable to reach its erosion target. But this second traverse almost always produced zero erosion, so seems not worth doing
//    int nRet = nTraversePolygonAndEstimateBeachErosion(nCoast, nPoly, dPotentialErosionOnPolygon, dTotFineEroded, dTotSandEroded, dTotCoarseEroded);
//    if (nRet != RTN_OK)
//       return nRet;
//
//    // Save these values
//    CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
//    pPolygon->SetDeltaEstimatedUnconsFine(-dTotFineEroded);
//    pPolygon->SetDeltaEstimatedUnconsSand(-dTotSandEroded);
//    pPolygon->SetDeltaEstimatedUnconsCoarse(-dTotCoarseEroded);
//
//    // Save the estimated values
//    m_dThisIterEstimatedActualFineBeachErosion   += dTotFineEroded;
//    m_dThisIterEstimatedActualSandBeachErosion   += dTotSandEroded;
//    m_dThisIterEstimatedActualCoarseBeachErosion += dTotCoarseEroded;
//
//    return RTN_OK;
// }

/*===============================================================================================================================

 This routine estimates how much we can erode on this polygon from profiles which are parallel to the polygon's up-coast boundary, moving down-coast (i.e. in the direction of increasing coastpoint indices). It does not do any actual erosion. The estimated values for fine, sand, and coarse erosion are depths on a raster cell: to convert to volumes, multiply by m_dCellArea

===============================================================================================================================*/
int CSimulation::nTraversePolygonAndEstimateBeachErosion(int const nCoast, int const nPoly, double const dErosionTarget, double &dTotFineEroded, double &dTotSandEroded, double &dTotCoarseEroded)
{
   double dStillToErodeOnPolygon = dErosionTarget;

   CGeomCoastPolygon *pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   // Get the up-coast and down-coast boundary details
   int nUpCoastProfile = pPolygon->nGetUpCoastProfile();
   CGeomProfile *pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);

   int nDownCoastProfile = pPolygon->nGetDownCoastProfile();
   CGeomProfile *pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We will use only part of the up-coast boundary profile, seaward as far as the depth of closure. First find the seaward end point of this up-coast part-profile. Note that this does not change as the landwards offset changes
   int nIndex = pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
   if (nIndex == INT_NODATA)
   {
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": " << ERR << "while estimating beach erosion for polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

      return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_1;
   }

   // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at the depth of closure
   int nUpCoastPartProfileLen = nIndex + 1;

   //    assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Store the cell co-ordinates of the boundary part-profile in reverse (sea to coast) order so we can append to the coastward end as we move inland (i.e. as nInlandOffset increases)
   vector<CGeom2DIPoint> PtiVUpCoastPartProfileCell;
   for (int n = 0; n < nUpCoastPartProfileLen; n++)
      PtiVUpCoastPartProfileCell.push_back(*pUpCoastProfile->pPtiGetCellInProfile(nUpCoastPartProfileLen - n - 1));

   int
       nUpCoastProfileCoastPoint = pUpCoastProfile->nGetNumCoastPoint(),
       nDownCoastProfileCoastPoint = pDownCoastProfile->nGetNumCoastPoint(),
       nXUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetX(),
       nYUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetY(),
       nCoastSegLen;

   // Store the coast point numbers for this polygon so that we can shuffle them
   vector<int> nVCoastPoint;
   if (nDownCoastProfileCoastPoint == m_VCoast[nCoast].nGetCoastlineSize() - 1)
   {
      // This is the final down-coast polygon, so also include the down-coast polygon boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }
   else
   {
      // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }

   // Estimate the volume of sediment which is to be eroded from each parallel profile
   double dAllSedimentTargetPerProfile = dErosionTarget / nCoastSegLen;

   // If there is zero (or near zero) sediment to erode, then just return
   if (bFPIsEqual(dAllSedimentTargetPerProfile, 0, TOLERANCE))
      return RTN_OK;

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

   // Traverse the polygon's existing coastline in a DOWN-COAST (i.e. increasing coastpoint indices) sequence, at each coast point fitting a Dean profile which is parallel to the up-coast polygon boundary
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Get the coast point
      int nCoastPoint = nVCoastPoint[n];

      CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int
          nCoastX = PtiCoastPoint.nGetX(),
          nCoastY = PtiCoastPoint.nGetY();

      //       LogStream << m_ulIter << ": nCoastX = " << nCoastX << ", nCoastY = " << nCoastY << ", this is {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}" << endl;

      // Is the coast cell an intervention structure?
      if (m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
      {
         // No erosion possible on this parallel profile, so move on
         //          LogStream << m_ulIter << ": intervention structure at coast point [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}, cannot erode this parallel profile" << endl;

         continue;
      }

      //       LogStream << m_ulIter << ": PtiVUpCoastPartProfileCell.back().nGetX() = " << PtiVUpCoastPartProfileCell.back().nGetX() << ", PtiVUpCoastPartProfileCell.back().nGetY() = " << PtiVUpCoastPartProfileCell.back().nGetY() << ", this is {" << dGridCentroidXToExtCRSX(PtiVUpCoastPartProfileCell.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVUpCoastPartProfileCell.back().nGetY()) << "}" << endl;

      // Not an intervention structure, so calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int
          nXOffset = nCoastX - PtiVUpCoastPartProfileCell.back().nGetX(),
          nYOffset = nCoastY - PtiVUpCoastPartProfileCell.back().nGetY();

      // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in reverse sequence, like the boundary part-profile)
      vector<CGeom2DIPoint> VPtiParProfile;
      for (int m = 0; m < nUpCoastPartProfileLen; m++)
      {
         // TODO check that each point is within valid grid, do same for other similar places in rest of model
         CGeom2DIPoint PtiTmp(PtiVUpCoastPartProfileCell[m].nGetX() + nXOffset, PtiVUpCoastPartProfileCell[m].nGetY() + nYOffset);
         VPtiParProfile.push_back(PtiTmp);
      }

      // Get the elevations of the start and end points of the parallel profiles (as we extend the profile inland, the elevation of the new coast point of the Dean profile is set to the elevation of the original coast point)
      int
          nParProfEndX = VPtiParProfile[0].nGetX(),
          nParProfEndY = VPtiParProfile[0].nGetY();

      // Safety check
      if (!bIsWithinValidGrid(nParProfEndX, nParProfEndY))
      {
         //          LogStream << WARN << "01 @@@@ while estimating actual beach erosion for coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
         VPtiParProfile[0].SetX(nParProfEndX);
         VPtiParProfile[0].SetY(nParProfEndY);
      }

      bool
          bHitEdge = false,
          bEndProfile = false,
          bZeroGradient = false,
          bEnoughEroded = false;

      int
          nParProfLen,
          nInlandOffset = -1;

      double
          dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev(),
          dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetSedimentTopElev();

      vector<double> VdParProfileDeanElev;

      // These are for saving values for each offset
      vector<int> VnParProfLenEachOffset;
      vector<double> VdAmountEachOffset;
      vector<vector<CGeom2DIPoint>> VVPtiParProfileEachOffset;
      vector<vector<double>> VVdParProfileDeanElevEachOffset;

      // OK, loop either until we can erode sufficient unconsolidated sediment, or until the landwards-moving parallel profile hits the grid edge
      while (true)
      {
         // Move inland by one cell
         nInlandOffset++;

         if (nInlandOffset > 0)
         {
            if (nInlandOffset > (pUpCoastProfile->nGetNumCellsInProfile() - 1))
            {
               //                LogStream << m_ulIter << ": reached end of up-coast profile " << nUpCoastProfile << " during down-coast ESTIMATION of actual beach erosion for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               bEndProfile = true;
               break;
            }

            // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing up-coast part-profile start point
            CGeom2DIPoint PtiUpCoastTmp = *pUpCoastProfile->pPtiGetCellInProfile(nInlandOffset);

            // Then get the offset between this PtiUpCoastTmp cell and the existing up-coast part-profile start point, and use the reverse of this offset to get the co-ords of the cell that extends the existing up-coast part-profile landwards
            int
                nXUpCoastStartOffset = PtiUpCoastTmp.nGetX() - nXUpCoastProfileExistingCoastPoint,
                nYUpCoastStartOffset = PtiUpCoastTmp.nGetY() - nYUpCoastProfileExistingCoastPoint,
                nXUpCoastThisStart = nCoastX - nXUpCoastStartOffset,
                nYUpCoastThisStart = nCoastY - nYUpCoastStartOffset;

            // Is the new landwards point within the raster grid?
            if (!bIsWithinValidGrid(nXUpCoastThisStart, nYUpCoastThisStart))
            {
               // It isn't
               //                LogStream << WARN << "reached edge of grid at [" << nXUpCoastThisStart << "][" << nYUpCoastThisStart << "] = {" << dGridCentroidXToExtCRSX(nXUpCoastThisStart) << ", " << dGridCentroidYToExtCRSY(nYUpCoastThisStart) << "} during DOWN-COAST estimation of beach erosion for coast " << nCoast << " polygon " << nPoly << ", nCoastPoint = " << nCoastPoint << ", this is [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}, nInlandOffset = " << nInlandOffset << ")" << endl << endl;

               // TODO Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
               bHitEdge = true;
               break;
            }

            CGeom2DIPoint PtiThisUpCoastStart(nXUpCoastThisStart, nYUpCoastThisStart);

            // Calculate the co-ords of a possible new landwards cell for the parallel profile
            int
                nXParNew = nXUpCoastThisStart + nXOffset,
                nYParNew = nYUpCoastThisStart + nYOffset;

            // Safety check
            if (!bIsWithinValidGrid(nXParNew, nYParNew))
            {
               //                LogStream << WARN << "02 @@@@ while estimating actual beach erosion on coast " << nCoast << " polygon " << nPoly << " (nInlandOffset = " << nInlandOffset << "), outside valid grid at [" << nXParNew << "][" << nYParNew << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end, is now [";

               KeepWithinValidGrid(nCoastX, nCoastY, nXParNew, nYParNew);

               //                LogStream << "[" << nXParNew << "][" << nYParNew << "] = {" << dGridCentroidXToExtCRSX(nXParNew) << ", " <<  dGridCentroidYToExtCRSY(nYParNew) << "}" << endl;

               // Is this cell already in the parallel profile?
               if ((VPtiParProfile.back().nGetX() != nXParNew) || (VPtiParProfile.back().nGetY() != nYParNew))
               {
                  // It isn't, so append it to the parallel profile
                  CGeom2DIPoint PtiTmp(nXParNew, nYParNew);
                  VPtiParProfile.push_back(PtiTmp);
               }
            }
            else
            {
               // No problem, so just append this to the parallel profile
               CGeom2DIPoint PtiTmp(nXParNew, nYParNew);
               VPtiParProfile.push_back(PtiTmp);
            }
         }

         nParProfLen = static_cast<int>(VPtiParProfile.size());

         if (nParProfLen < MIN_PAR_PROFILE_SIZE)
         {
            // Can't have a meaningful parallel profile with very few points
            // LogStream << m_ulIter << ": only " << nParProfLen << " points in parallel profile, min is " << MIN_PAR_PROFILE_SIZE << ", abandoning" << endl;

            continue;
         }

         //          for (int m = 0; m < static_cast<int>(VPtiParProfile.size()); m++)
         //             LogStream << "[" << VPtiParProfile[m].nGetX() << "][" << VPtiParProfile[m].nGetY() << "] ";
         //          LogStream << endl;

         // Get the distance between the start and end of the parallel profile, in external CRS units. Note that the parallel profile co-ords are in reverse sequence
         CGeom2DPoint
             PtStart = PtGridCentroidToExt(&VPtiParProfile.back()),
             PtEnd = PtGridCentroidToExt(&VPtiParProfile[0]);

         // Calculate the length of the parallel profile
         double dParProfileLen = dGetDistanceBetween(&PtStart, &PtEnd);

         // Calculate the elevation difference between the start and end of the parallel profile
         double dElevDiff = dParProfCoastElev - dParProfEndElev;
         if (bFPIsEqual(dElevDiff, 0, TOLERANCE))
         {
            // Can't have a meaningful Dean profile with a near-zero elevation difference
            // TODO Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
            // LogStream << m_ulIter << ": zero gradient on parallel profile, abandoning" << endl;

            bZeroGradient = true;
            break;
         }

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double dParProfA = dElevDiff / pow(dParProfileLen, DEAN_POWER);
         VdParProfileDeanElev.resize(nParProfLen, 0);
         double dInc = dParProfileLen / (nParProfLen - 1);

         // For the parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfCoastElev, dParProfA, false, 0, 0);

         vector<double> dVParProfileNow(nParProfLen, 0);
         vector<bool> bVProfileValid(nParProfLen, true);
         for (int m = 0; m < nParProfLen; m++)
         {
            int
                nX = VPtiParProfile[nParProfLen - m - 1].nGetX(),
                nY = VPtiParProfile[nParProfLen - m - 1].nGetY();

            // Safety check
            if (!bIsWithinValidGrid(nX, nY))
            {
               //                LogStream << WARN << "03 @@@@ while constructing parallel profile to estimate DOWN-COAST actual beach erosion on coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               bVProfileValid[m] = false;
               continue;
            }

            // Don't erode intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               bVProfileValid[m] = false;

            dVParProfileNow[m] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();
         }

         // Get the total difference in elevation (present profile - Dean profile)
         double dParProfTotDiff = dSubtractProfiles(&dVParProfileNow, &VdParProfileDeanElev, &bVProfileValid);

         // DEBUG STUFF -----------------------------------------------------
         //          LogStream << m_ulIter<< ": in polygon " << nPoly << " at coast point " << nCoastPoint << " doing DOWN-COAST estimation of actual beach erosion, parallel profile with nInlandOffset = " << nInlandOffset << ", from [" << VPtiParProfile.back().nGetX() << "][" << VPtiParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile.back().nGetY()) << "} to [" << VPtiParProfile[0].nGetX() << "][" << VPtiParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;

         //          LogStream << "Profile now:" << endl;
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             if (bVProfileValid[n])
         //                LogStream << dVParProfileNow[n] << " ";
         //             else
         //                LogStream << "XXX ";
         //          }
         //          LogStream << endl << endl;;
         //
         //          LogStream << "Parallel Dean profile for erosion:" << endl;
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             if (bVProfileValid[n])
         //                LogStream << VdParProfileDeanElev[n] << " ";
         //             else
         //                LogStream << "XXX ";
         //          }
         //          LogStream << endl << endl;
         //
         //          LogStream << "Difference (present profile minus Dean profile):" << endl;
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             if (bVProfileValid[n])
         //                LogStream << dVParProfileNow[n] - VdParProfileDeanElev[n] << " ";
         //             else
         //                LogStream << "XXX ";
         //          }
         //          LogStream << endl << endl;
         //
         //          LogStream << "dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
         // DEBUG STUFF -----------------------------------------------------

         // So will we be able to erode as much as is needed?
         if (dParProfTotDiff > dAllSedimentTargetPerProfile)
         {
            //             LogStream << m_ulIter << ": in polygon " << nPoly << " at coast point " << nCoastPoint << " doing DOWN-COAST estimation of actual beach erosion for parallel profile with nInlandOffset = " << nInlandOffset << ", can meet erosion target: dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

            bEnoughEroded = true;
            break;
         }

         // We have not been able to reach the erosion target for this parallel profile. Now, if we have moved at least MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION cells inland, and dParProfTotDiff is zero, break out of the loop
         if ((nInlandOffset >= MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION) && (bFPIsEqual(dParProfTotDiff, 0, TOLERANCE)))
            break;

         // Save the amount which can be eroded for this offset TODO save other stuff too
         VdAmountEachOffset.push_back(dParProfTotDiff);
         VnParProfLenEachOffset.push_back(nParProfLen);
         VVPtiParProfileEachOffset.push_back(VPtiParProfile);
         VVdParProfileDeanElevEachOffset.push_back(VdParProfileDeanElev);
      }

      // If we hit the edge of the grid, or have a zero gradient on the profile, or this is an end profile, then abandon this profile and do the next parallel profile
      // TODO Improve this, see above
      if (bHitEdge || bEndProfile || bZeroGradient)
         continue;

      // OK we will do some erosion on this parallel profile, set the target for this profile to be the full target amount
      double dStillToErodeOnProfile = dAllSedimentTargetPerProfile;

      if (!bEnoughEroded)
      {
         // We have not been able to reach the target for erosion on this parallel profile. So find the offset that gives us the largest erosion amount
         int nOffsetForLargestPossible = -1;
         double dLargestPossibleErosion = 0;
         for (unsigned int nn = 0; nn < VdAmountEachOffset.size(); nn++)
         {
            if (VdAmountEachOffset[nn] > dLargestPossibleErosion)
            {
               dLargestPossibleErosion = VdAmountEachOffset[nn];
               nOffsetForLargestPossible = nn;
            }
         }

         // If every offset gave zero erosion then abandon this profile and do the next parallel profile
         if (nOffsetForLargestPossible < 0)
            continue;

         // OK, we have an offset which gives us the largest possible erosion (but less than the full target amount), continue with this
         nInlandOffset = nOffsetForLargestPossible;
         dStillToErodeOnProfile = dLargestPossibleErosion;
         nParProfLen = VnParProfLenEachOffset[nInlandOffset];
         VPtiParProfile = VVPtiParProfileEachOffset[nInlandOffset];
         VdParProfileDeanElev = VVdParProfileDeanElevEachOffset[nInlandOffset];

         //          LogStream << m_ulIter << ": in polygon " << nPoly << " at coast point " << nCoastPoint << " doing DOWN-COAST estimation of actual beach erosion for parallel profile with nInlandOffset = " << nInlandOffset << ", could not meet erosion target dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << ", instead using best possible: nInlandOffset = " << nInlandOffset << " which gives dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;
      }

      // OK, this value of nInlandOffset gives us some (tho' maybe not enough) erosion. So calculate how much fine-, sand-, and coarse sized sediment can be obtained here by working along the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0)
      int nRet = nEstimateBeachErosionOnParallelProfile(/*nPoly, nCoastPoint,*/ nCoastX, nCoastY, /* nInlandOffset, */ nParProfLen, &VPtiParProfile, &VdParProfileDeanElev, dStillToErodeOnProfile, dStillToErodeOnPolygon, dTotFineEroded, dTotSandEroded, dTotCoarseEroded);
      if (nRet != RTN_OK)
         return nRet;

      //       LogStream << m_ulIter << ": in polygon " << nPoly << " estimating erosion: finished at coast point " << nCoastPoint << " dStillToErodeOnProfile = " << dStillToErodeOnProfile << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " amount eroded = " << dAllSedimentTargetPerProfile - dStillToErodeOnProfile << endl;
   }

   return RTN_OK;
}

/*===============================================================================================================================

 This routine calculates the supply-limited amount of beach erosion, in three size categories, which can be eroded on this parallel profile. It does not do any actual erosion

 ===============================================================================================================================*/
int CSimulation::nEstimateBeachErosionOnParallelProfile(/*int const nPoly, int const nCoastPoint,*/ int const nCoastX, int const nCoastY, /* int const nInlandOffset, */ int const nParProfLen, vector<CGeom2DIPoint> const *pVPtiParProfile, vector<double> const *pVdParProfileDeanElev, double &dStillToErodeOnProfile, double &dStillToErodeOnPolygon, double &dTotFineEroded, double &dTotSandEroded, double &dTotCoarseEroded)
{
   double
       dFineEroded = 0, // Totals for this parallel profile
       dSandEroded = 0,
       dCoarseEroded = 0;

   for (int nDistSeawardFromNewCoast = 0; nDistSeawardFromNewCoast < nParProfLen; nDistSeawardFromNewCoast++)
   {
      // Don't bother with tiny amounts
      if (dStillToErodeOnProfile < SEDIMENT_ELEV_TOLERANCE)
         dStillToErodeOnProfile = 0;

      // Leave the loop if we have eroded enough for this polygon
      if (dStillToErodeOnPolygon <= 0)
      {
         //          LogStream << m_ulIter<< ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for polygon, dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << endl;

         break;
      }

      // Leave the loop if we have done enough erosion for this profile
      if (dStillToErodeOnProfile <= 0)
      {
         //          LogStream << m_ulIter << ": in polygon " << nPoly << " at coast point " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for profile, dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;

         break;
      }

      CGeom2DIPoint PtiTmp = pVPtiParProfile->at(nParProfLen - nDistSeawardFromNewCoast - 1);
      int
          nX = PtiTmp.nGetX(),
          nY = PtiTmp.nGetY();

      // Safety check
      if (!bIsWithinValidGrid(nX, nY))
      {
         //             LogStream << WARN << "04 @@@@ while estimating sediment size fractions on eroding polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
         PtiTmp.SetX(nX);
         PtiTmp.SetY(nY);
      }

      // Don't do anything to intervention cells
      if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
         continue;

      // Don't do cells twice
      if (!m_pRasterGrid->m_Cell[nX][nY].bGetActualBeachErosionEstimated())
      {
         // Set flag so we don't estimate erosion/deposition on this cell again
         m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosionEstimated();

         // Get this cell's current elevation
         double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

         //             LogStream << "\tnPoly = " << nPoly << ", [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

         // Subtract the two elevations
         double dElevDiff = dThisElevNow - pVdParProfileDeanElev->at(nDistSeawardFromNewCoast);
         if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
         {
            // The current elevation is higher than the Dean elevation, so we have possible beach erosion (i.e. if not constrained by availability of unconsolidated sediment) here
            //                LogStream << "\tnPoly = " << nPoly << " doing DOWN-COAST, possible beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

            // Now get the number of the highest layer with non-zero thickness
            int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

            // Safety check
            if (nThisLayer == INT_NODATA)
               return RTN_ERR_NO_TOP_LAYER;

            if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
            {
               // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So estimate some beach erosion
               double
                   dToErode = tMin(dElevDiff, dStillToErodeOnProfile, dStillToErodeOnPolygon),
                   dFine = 0,
                   dSand = 0,
                   dCoarse = 0;

               //                assert(dToErode > 0);

               EstimateActualBeachErosionOnCell(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

               // Totals for this parallel profile
               dFineEroded += dFine;
               dSandEroded += dSand;
               dCoarseEroded += dCoarse;

               dTotFineEroded += dFine;
               dTotSandEroded += dSand;
               dTotCoarseEroded += dCoarse;

               double dTmpTot = dFine + dSand + dCoarse;
               dStillToErodeOnProfile -= dTmpTot;
               dStillToErodeOnPolygon -= dTmpTot;

               //                   LogStream << m_ulIter << ": in polygon" << nPoly << ", estimated actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
            }
         }
         else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
         {
            // The current elevation is below the Dean elevation, so we have can have beach deposition here provided that we have some previously-eroded unconsolidated sediment (sand and coarse only) to deposit
            double dTotSandAndCoarseEroded = dTotSandEroded + dTotCoarseEroded;
            if (dTotSandAndCoarseEroded > SEDIMENT_ELEV_TOLERANCE)
            {
               // Assume that the sediment size fractions which are deposited are in the same ratios as the previously-eroded sediment
               double
                   dToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                   dSandProp = dTotSandEroded / dTotSandAndCoarseEroded,
                   dCoarseProp = 1 - dSandProp,
                   dSandToDeposit = tMin(dToDeposit * dSandProp, dTotSandEroded),
                   dCoarseToDeposit = tMin(dToDeposit * dCoarseProp, dTotCoarseEroded);

               if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
               {
                  // Totals for this parallel profile
                  dSandEroded -= dSandToDeposit;
                  dStillToErodeOnProfile += dSandToDeposit;

                  // Totals for the polygon
                  dTotSandEroded -= dSandToDeposit;
                  dStillToErodeOnPolygon += dSandToDeposit;
               }

               if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
               {
                  // Totals for this parallel profile
                  dCoarseEroded -= dCoarseToDeposit;
                  dStillToErodeOnProfile += dCoarseToDeposit;

                  // Totals for the polygon
                  dTotCoarseEroded -= dCoarseToDeposit;
                  dStillToErodeOnPolygon += dCoarseToDeposit;
               }

               //                LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << ", beach deposition = " << dSandToDeposit + dCoarseToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}" << endl;
            }
         }
      }
   }

   return RTN_OK;
}

/*===============================================================================================================================

 Estimates actual (supply-limited) erosion of unconsolidated beach sediment on a single cell, returns the depth-equivalents of fine, sand and coarse sediment which could be removed

===============================================================================================================================*/
void CSimulation::EstimateActualBeachErosionOnCell(int const nX, int const nY, int const nThisLayer, double const dPotentialErosionOnPolygon, double &dFine, double &dSand, double &dCoarse)
{
   // Find out how much unconsolidated sediment we have available on this cell
   double
       dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFine(),
       dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSand(),
       dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarse();

   // Is there any unconsolidated sediment on this cell?
   if ((dExistingAvailableFine + dExistingAvailableSand + dExistingAvailableCoarse) < SEDIMENT_ELEV_TOLERANCE)
      return;

   // We have some unconsolidated sediment, so partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   int
       nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0),
       nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0),
       nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

   double dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);

   if (nFineWeight)
   {
      // Pretend to erode some fine-sized consolidated sediment
      double dFineLowering = (m_dFineErodibilityNormalized * dPotentialErosionOnPolygon) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      dFine = tMin(dExistingAvailableFine, dFineLowering);
   }

   if (nSandWeight)
   {
      // Pretend to erode some sand-sized consolidated sediment
      double dSandLowering = (m_dSandErodibilityNormalized * dPotentialErosionOnPolygon) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSand = tMin(dExistingAvailableSand, dSandLowering);
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double dCoarseLowering = (m_dCoarseErodibilityNormalized * dPotentialErosionOnPolygon) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarse = tMin(dExistingAvailableCoarse, dCoarseLowering);
   }

   //    LogStream << "\tEstimating actual erosion on [" << nX << "][" << nY << "] dPotentialErosionOnPolygon = " << dPotentialErosionOnPolygon << " dTotActualErosion = " << dTotActualErosion << " dFine = " << dFine << " dSand = " << dSand << " dCoarse = " << dCoarse << endl;
}
