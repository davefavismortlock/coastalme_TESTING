/*!
 *
 * \file do_beach_within_polygon.cpp
 * \brief Does within-polygon actual erosion and distribution of transported beach sediment
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
// #include <assert.h>

#include <cmath>
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

/*===============================================================================================================================

 Erodes unconsolidated beach sediment on the cells within a polygon

 As with the estimation of actual erosion, this is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For erosion, this reduces the beach volume

===============================================================================================================================*/
int CSimulation::nDoBeachErosionOnPolygon(int const nCoast, int const nPoly, double const dErosionTarget, double &dFineErosionDiff, double &dSandErosionDiff, double &dCoarseErosionDiff)
{
   double
       dTotFineEroded = 0,
       dTotSandEroded = 0,
       dTotCoarseEroded = 0;

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
      LogStream << m_ulIter << ": " << ERR << "while eroding beach in polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

      return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_2;
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
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << WARN << "while eroding beach erosion for coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

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
               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << m_ulIter << ": reached end of up-coast profile " << nUpCoastProfile << " during down-coast beach erosion for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

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
               //                LogStream << WARN << "reached edge of grid at [" << nXUpCoastThisStart << "][" << nYUpCoastThisStart << "] = {" << dGridCentroidXToExtCRSX(nXUpCoastThisStart) << ", " << dGridCentroidYToExtCRSY(nYUpCoastThisStart) << "} during DOWN-COAST beach erosion for coast " << nCoast << " polygon " << nPoly << ", nCoastPoint = " << nCoastPoint << ", this is [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " <<  dGridCentroidYToExtCRSY(nCoastY) << "}, nInlandOffset = " << nInlandOffset << ")" << endl << endl;

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
               //                LogStream << WARN << "while eroding beach on coast " << nCoast << " polygon " << nPoly << " (nInlandOffset = " << nInlandOffset << "), outside valid grid at [" << nXParNew << "][" << nYParNew << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end, is now [";

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
               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << WARN << "while constructing parallel profile for beach erosion on coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

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
         //          LogStream << m_ulIter<< ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", parallel profile with nInlandOffset = " << nInlandOffset << ", from [" << VPtiParProfile.back().nGetX() << "][" << VPtiParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile.back().nGetY()) << "} to [" << VPtiParProfile[0].nGetX() << "][" << VPtiParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(VPtiParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;

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
            //             LogStream << m_ulIter << ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", on parallel profile with nInlandOffset = " << nInlandOffset << ", can meet erosion target: dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

            bEnoughEroded = true;
            break;
         }

         // We have not been able to reach the erosion target for this parallel profile. Now, if we have moved at least MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION cells inland, and dParProfTotDiff is zero, break out of the loop
         if ((nInlandOffset >= MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION) && (bFPIsEqual(dParProfTotDiff, 0, TOLERANCE)))
         {
            if (m_nLogFileDetail >= LOG_FILE_MOST_DETAIL)
               LogStream << m_ulIter << ": leaving loop because nInlandOffset (" << nInlandOffset << ") >= MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION) and dParProfTotDiff = " << dParProfTotDiff << endl;
            break;
         }

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

         //          LogStream << m_ulIter << ": eroding polygon " << nPoly << " at coast point " << nCoastPoint << ", for parallel profile with nInlandOffset = " << nInlandOffset << ", could not meet erosion target dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << ", instead using best possible: nInlandOffset = " << nInlandOffset << " which gives dStillToErodeOnProfile = " << dStillToErodeOnProfile << endl;
      }

      // This value of nInlandOffset gives us some (tho' maybe not enough) erosion. So do the erosion, by working along the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0)
      int nRet = nDoBeachErosionOnParallelProfile(/* nPoly, nCoastPoint, */ nCoastX, nCoastY, /* nInlandOffset, */ nParProfLen, &VPtiParProfile, &VdParProfileDeanElev, dStillToErodeOnProfile, dStillToErodeOnPolygon, dTotFineEroded, dTotSandEroded, dTotCoarseEroded);
      if (nRet != RTN_OK)
         return nRet;

      //       LogStream << m_ulIter << ": eroding polygon " << nPoly << ": finished at coast point " << nCoastPoint << " dStillToErodeOnProfile = " << dStillToErodeOnProfile << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

      if (dStillToErodeOnProfile > 0)
      {
         //          LogStream << "XXXXXXXXXXXXXXXXXXXXXX Still to erode on profile = " << dStillToErodeOnProfile << " dStillToErodeOnPolygon WAS = " << dStillToErodeOnPolygon;
         dStillToErodeOnPolygon += dStillToErodeOnProfile;
         dAllSedimentTargetPerProfile += (dStillToErodeOnProfile / (nCoastSegLen - n));
         //          LogStream << " dStillToErodeOnPolygon NOW = " << dStillToErodeOnPolygon << endl;
      }

      if (dStillToErodeOnPolygon <= 0)
      {
         //          LogStream << m_ulIter << ": YYYYYYYYYYYYYYY in polygon " << nPoly << ", leaving loop because dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << endl;
         break;
      }
   }

   // How much have we been able to erode?
   //    LogStream << m_ulIter << ": nPoly = " << nPoly << " dTotFineEroded = " << dTotFineEroded << " dTotSandEroded = " << dTotSandEroded << " dTotCoarseEroded = " << dTotCoarseEroded << endl;

   double dTotalErosion = dTotFineEroded + dTotSandEroded + dTotCoarseEroded;

   // Is the depth eroded within TOLERANCE of the target depth-equivalent?
   if (bFPIsEqual(dTotalErosion, dErosionTarget, TOLERANCE))
   {
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": polygon " << nPoly << " actual beach erosion approx equal to target beach erosion: actual = " << dTotalErosion << " target = " << dErosionTarget << endl;
   }
   else
   {
      if (dTotalErosion < dErosionTarget)
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": on polygon " << nPoly << " actual beach erosion is less than target beach erosion: actual = " << dTotalErosion << " target = " << dErosionTarget << " difference = " << -(dErosionTarget - dTotalErosion) << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << ", will reduce sediment exported from this polygon" << endl;
      }
      else
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << ERR << "on polygon " << nPoly << " actual beach erosion is GREATER THAN target beach erosion: actual = " << dTotalErosion << " target = " << dErosionTarget << " difference = " << -(dErosionTarget - dTotalErosion) << " dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << endl;
      }

      // Calculate the difference between what we should have been able to erode, and what we actually eroded
      double dTotErosionDiff = dErosionTarget - dTotalErosion;
      // LogStream << endl << "&&&0 nCoast = " << nCoast << " nPoly = " << nPoly << " dTotErosionDiff = " << dTotErosionDiff << endl;

      if (dTotalErosion > 0)
      {
         // Distribute this difference (a +ve number) between potential and actual between the size classes
         dFineErosionDiff = dTotErosionDiff * dTotFineEroded / dTotalErosion;
         dSandErosionDiff = dTotErosionDiff * dTotSandEroded / dTotalErosion;
         dCoarseErosionDiff = dTotErosionDiff * dTotCoarseEroded / dTotalErosion;

         m_dThisIterErosionFineDiff += dFineErosionDiff;
         m_dThisIterErosionSandDiff += dSandErosionDiff;
         m_dThisIterErosionCoarseDiff += dCoarseErosionDiff;

         // LogStream << "&&&1 dFineErosionDiff = " << dFineErosionDiff << " m_dThisIterErosionFineDiff = " << m_dThisIterErosionFineDiff << endl;
         // LogStream << "&&&1 dSandErosionDiff = " << dSandErosionDiff << " m_dThisIterErosionSandDiff = " << m_dThisIterErosionSandDiff  << endl;
         // LogStream << "&&&1 dCoarseErosionDiff = " << dCoarseErosionDiff << " m_dThisIterErosionCoarseDiff = " << m_dThisIterErosionCoarseDiff  << endl;
      }
   }

   return RTN_OK;
}

/*===============================================================================================================================

 This routine erodes the unconsolidated beach sediment on this parallel profile

===============================================================================================================================*/
int CSimulation::nDoBeachErosionOnParallelProfile(/* int const nPoly, int const nCoastPoint, */ int const nCoastX, int const nCoastY, /* int const nInlandOffset, */ int const nParProfLen, vector<CGeom2DIPoint> const *pVPtiParProfile, vector<double> const *pVdParProfileDeanElev, double &dStillToErodeOnProfile, double &dStillToErodeOnPolygon, double &dTotFineEroded, double &dTotSandEroded, double &dTotCoarseEroded)
{
   double
       dFineEroded = 0, // Totals for this parallel profile
       dSandEroded = 0,
       dCoarseEroded = 0;

   for (int nDistSeawardFromNewCoast = 0; nDistSeawardFromNewCoast < nParProfLen; nDistSeawardFromNewCoast++)
   {
      // Leave the loop if we have eroded enough for this polygon
      if (dStillToErodeOnPolygon <= 0)
      {
         //          LogStream << m_ulIter<< ": in polygon " << nPoly << " at coast point " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for polygon, dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << " (dStillToErodeOnProfile = " << dStillToErodeOnProfile << ")" << endl;

         break;
      }

      // Leave the loop if we have done enough erosion for this profile
      if (dStillToErodeOnProfile <= 0)
      {
         //          LogStream << m_ulIter << ": in polygon " << nPoly << " at coast point " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ", leaving loop because enough erosion for profile, dStillToErodeOnProfile = " << dStillToErodeOnProfile << " (dStillToErodeOnPolygon = " << dStillToErodeOnPolygon << ")" << endl;

         break;
      }

      CGeom2DIPoint PtiTmp = pVPtiParProfile->at(nParProfLen - nDistSeawardFromNewCoast - 1);
      int
          nX = PtiTmp.nGetX(),
          nY = PtiTmp.nGetY();

      // Safety check
      if (!bIsWithinValidGrid(nX, nY))
      {
         //          LogStream << WARN << "04 @@@@ while eroding polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << ". Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
         PtiTmp.SetX(nX);
         PtiTmp.SetY(nY);
      }

      // Don't do anything to intervention cells
      if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
         continue;

      // Don't do cells twice
      if (!m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
      {
         // Get this cell's current elevation
         double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

         //             LogStream << "\tnPoly = " << nPoly << ", [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

         // Subtract the two elevations
         double dElevDiff = dThisElevNow - pVdParProfileDeanElev->at(nDistSeawardFromNewCoast);
         if ((dElevDiff > 0) && (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())) // added the check of contiguos sea MCB
         {
            // The current elevation is higher than the Dean elevation, so we have possible beach erosion (i.e. if not constrained by availability of unconsolidated sediment) here
            //             LogStream << "\tnPoly = " << nPoly << " doing DOWN-COAST, possible beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << pVdParProfileDeanElev->at(nDistSeawardFromNewCoast) << endl;

            // Now get the number of the highest layer with non-zero thickness
            int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

            // Safety check
            if (nThisLayer == INT_NODATA)
               return RTN_ERR_NO_TOP_LAYER;

            if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
            {
               // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
               double
                   dToErode = tMin(dElevDiff, dStillToErodeOnProfile, dStillToErodeOnPolygon),
                   dFine = 0,
                   dSand = 0,
                   dCoarse = 0;

               //                assert(dToErode > 0);

               ErodeBeachSedimentOnCellSupplyLimited(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

               double dTmpTot = dFine + dSand + dCoarse;
               if (dTmpTot > 0)
               {
                  // Update totals for this parallel profile
                  dFineEroded += dFine;
                  dSandEroded += dSand;
                  dCoarseEroded += dCoarse;

                  // Update totals for the polygon
                  dTotFineEroded += dFine;
                  dTotSandEroded += dSand;
                  dTotCoarseEroded += dCoarse;

                  dStillToErodeOnProfile -= dTmpTot;
                  dStillToErodeOnPolygon -= dTmpTot;

                  // Update this-timestep totals
                  m_ulThisIterNumActualBeachErosionCells++;
                  m_dThisIterActualBeachErosionFine += dFine;
                  m_dThisIterActualBeachErosionSand += dSand;
                  m_dThisIterActualBeachErosionCoarse += dCoarse;

                  //                   LogStream << m_ulIter << ": in polygon " << nPoly << ", actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
         }
         else if ((dElevDiff < 0) && (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()))
         {
            // The current elevation is below the Dean elevation, so we have can have beach deposition here provided that we have some previously-eroded unconsolidated sediment (sand and coarse only) to deposit
            double dTotSandAndCoarseEroded = dTotSandEroded + dTotCoarseEroded;
            if (dTotSandAndCoarseEroded > SEDIMENT_ELEV_TOLERANCE)
            {
               // Assume that the sediment size fractions which are deposited are in the same ratios as the previously-eroded sediment
               double
                   dTotToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                   dSandToDeposit = dTotToDeposit * dTotSandEroded / dTotSandAndCoarseEroded,
                   dCoarseToDeposit = dTotToDeposit * dTotCoarseEroded / dTotSandAndCoarseEroded;

               int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

               // Safety check
               if (nTopLayer == INT_NODATA)
                  return RTN_ERR_NO_TOP_LAYER;

               if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
               {
                  dSandToDeposit = tMin(dSandToDeposit, dTotSandEroded);

                  double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotSandEroded -= dSandToDeposit;

                  dStillToErodeOnProfile += dSandToDeposit;
                  dStillToErodeOnPolygon += dSandToDeposit;
               }

               if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
               {
                  dCoarseToDeposit = tMin(dCoarseToDeposit, dTotCoarseEroded);

                  double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotCoarseEroded -= dCoarseToDeposit;

                  dStillToErodeOnProfile += dCoarseToDeposit;
                  dStillToErodeOnPolygon += dCoarseToDeposit;
               }

               // Now update the cell's layer elevations
               m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

               // Update the cell's sea depth
               m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

               // Update the cell's beach deposition, and total beach deposition, values
               m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit + dCoarseToDeposit);

               // And set the landform category
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

               // Update per-timestep totals
               m_ulThisIterNumBeachDepositionCells++;
               m_dThisIterBeachDepositionSand += dSandToDeposit;
               m_dThisIterBeachDepositionCoarse += dCoarseToDeposit;

               //             LogStream << m_ulIter << ": nPoly = " << nPoly << ", beach deposition = " << dTotToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
            }
         }
      }
   }

   return RTN_OK;
}

/*===============================================================================================================================

 Erodes the unconsolidated beach sediment on a single cell, returns the depth-equivalents of fine, sand and coarse sediment removed

===============================================================================================================================*/
void CSimulation::ErodeBeachSedimentOnCellSupplyLimited(int const nX, int const nY, int const nThisLayer, double const dMaxToErode, double &dFine, double &dSand, double &dCoarse)
{
   // Find out how much unconsolidated sediment we have available on this cell
   double
       dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFine(),
       dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSand(),
       dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarse();

   // Is there any unconsolidated sediment on this cell?
   if ((dExistingAvailableFine + dExistingAvailableSand + dExistingAvailableCoarse) <= 0)
   {
      return;
   }

   // We have some unconsolidated sediment, so partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   int
       nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0),
       nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0),
       nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

   double
       dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized),
       dTotActualErosion = 0;

   if (nFineWeight)
   {
      // Erode some fine-sized consolidated sediment
      double dFineLowering = (m_dFineErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      dFine = tMin(dExistingAvailableFine, dFineLowering);
      double dRemaining = dExistingAvailableFine - dFine;

      dTotActualErosion += dFine;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetFine(dRemaining);

      // And set the changed-this-timestep switch
      m_bUnconsChangedThisIter[nThisLayer] = true;
   }

   if (nSandWeight)
   {
      // Erode some sand-sized consolidated sediment
      double dSandLowering = (m_dSandErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSand = tMin(dExistingAvailableSand, dSandLowering);
      double dRemaining = dExistingAvailableSand - dSand;

      dTotActualErosion += dSand;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetSand(dRemaining);

      // Set the changed-this-timestep switch
      m_bUnconsChangedThisIter[nThisLayer] = true;
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double dCoarseLowering = (m_dCoarseErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarse = tMin(dExistingAvailableCoarse, dCoarseLowering);
      double dRemaining = dExistingAvailableCoarse - dCoarse;

      dTotActualErosion += dCoarse;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetCoarse(dRemaining);

      // Set the changed-this-timestep switch
      m_bUnconsChangedThisIter[nThisLayer] = true;
   }

   // Set the actual erosion value for this cell
   m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosion(dTotActualErosion);

   if (dTotActualErosion > 0)
   {
      // Recalculate the elevation of every layer
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

      // And update the cell's sea depth
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
   }

   //    LogStream << "\t[" << nX << "][" << nY << "] dPotentialErosion = " << dPotentialErosion << " dTotActualErosion = " << dTotActualErosion << " dFine = " << dFine << " dSand = " << dSand << " dCoarse = " << dCoarse << endl;
}

/*===============================================================================================================================

 Deposits unconsolidated beach sediment on the cells within a polygon

 As with the estimation of actual erosion, this is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For deposition, this adds to the beach volume

===============================================================================================================================*/
int CSimulation::nDoBeachDepositionOnPolygon(int const nCoast, int const nPoly, double const dTotTargetToDeposit)
{
   // Don't bother with tiny amounts
   if (dTotTargetToDeposit < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   CGeomCoastPolygon *pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   // Get the grid cell co-ords of this polygon's up-coast and down-coast profiles
   int
       nUpCoastProfile = pPolygon->nGetUpCoastProfile(),
       nDownCoastProfile = pPolygon->nGetDownCoastProfile();

   CGeomProfile *pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);
   CGeomProfile *pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We are using only part of each profile, seaward as far as the depth of closure. First find the seaward end point of the up-coast part-profile
   //    CGeom2DIPoint PtiUpCoastPartProfileSeawardEnd;
   //    int nIndex =  pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiUpCoastPartProfileSeawardEnd);
   int nIndex = pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
   if (nIndex == INT_NODATA)
   {
      LogStream << m_ulIter << ": " << ERR << "while depositing beach for coast " << nCoast << " polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;

      return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_3;
   }

   // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
   int nUpCoastDeanLen = nIndex + 1;

   //    assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
   //    CGeom2DPoint
   //       PtUpCoastProfileStart = *pUpCoastProfile->pPtGetPointInProfile(0),
   //       PtUpCoastProfileEnd = PtGridCentroidToExt(&PtiUpCoastPartProfileSeawardEnd);

   //   double dUpCoastDeanLen = dGetDistanceBetween(&PtUpCoastProfileStart, &PtUpCoastProfileEnd);

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

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

   // Get the volume of sediment which is to be deposited on the polygon and on each parallel profile
   double
       dSandToDepositOnPoly = pPolygon->dGetDeltaActualUnconsSand(),
       dCoarseToDepositOnPoly = pPolygon->dGetDeltaActualUnconsCoarse(),
       dAllSedimentTargetPerProfile = dTotTargetToDeposit / nCoastSegLen,
       dSandTargetPerProfile = dSandToDepositOnPoly / nCoastSegLen,
       dCoarseTargetPerProfile = dCoarseToDepositOnPoly / nCoastSegLen,
       dTotSandDeposited = 0, // Grand totals for the whole polygon
       dTotCoarseDeposited = 0;

   // Now traverse the polygon's existing coastline in a random (but broadly DOWN-COAST i.e. increasing coast point indices) sequence, fitting a Dean profile at each coast point
   // DOWN-COAST (increasing coastpoint indices) ================================================================================
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Pick a random coast point
      int nCoastPoint = nVCoastPoint[n];
      CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int
          nCoastX = PtiCoastPoint.nGetX(),
          nCoastY = PtiCoastPoint.nGetY();

      // Calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int
          nXOffset = nCoastX - nXUpCoastProfileExistingCoastPoint,
          nYOffset = nCoastY - nYUpCoastProfileExistingCoastPoint,
          nSeawardOffset = -1;
      //          nParProfLen;
      vector<CGeom2DIPoint> PtiVParProfile;
      vector<double> VdParProfileDeanElev;

      // OK, loop until we can deposit sufficient unconsolidated sediment
      while (true)
      {
         // Move seaward by one cell
         nSeawardOffset++;

         // And lengthen the parallel profile
         int nParProfLen = nUpCoastDeanLen + nSeawardOffset;

         if (nParProfLen > (pUpCoastProfile->nGetNumCellsInProfile()))
         {
            // We've reached the seaward end of the up-coast profile, need to quit
            // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            //    LogStream << m_ulIter << ": " << WARN << "reached seaward end of up-coast profile during DOWN-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

            break;
         }

         // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
         PtiVParProfile.resize(0);
         for (int m = 0; m < nParProfLen; m++)
         {
            CGeom2DIPoint PtiProf = *pUpCoastProfile->pPtiGetCellInProfile(m);
            CGeom2DIPoint PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
            PtiVParProfile.push_back(PtiTmp);
         }

         // Get the existing elevation of the seaward end of the parallel profile
         int
             nSeaEndX = PtiVParProfile.back().nGetX(),
             nSeaEndY = PtiVParProfile.back().nGetY();

         // Safety check
         if (!bIsWithinValidGrid(nSeaEndX, nSeaEndY))
         {
            if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               LogStream << WARN << "09 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
            PtiVParProfile.back().SetX(nSeaEndX);
            PtiVParProfile.back().SetY(nSeaEndY);
         }

         double dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetSedimentTopElev();

         // Set the start elevation for the Dean profile just a bit above SWL for this timestep, so that it is a Bruun profile
         double dParProfStartElev = m_dThisIterSWL + m_dDeanProfileStartAboveSWL;

         // Calculate the total length of the parallel profile, including any seaward offset
         double dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

         // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
         double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double dParProfA = (dParProfStartElev - dParProfEndElev) / pow(dParProfDeanLen, DEAN_POWER);

         nParProfLen = static_cast<int>(PtiVParProfile.size());
         VdParProfileDeanElev.resize(nParProfLen, 0);

         //          for (int m = 0; m < static_cast<int>(PtiVParProfile.size()); m++)
         //             LogStream << "[" << PtiVParProfile[m].nGetX() << "][" << PtiVParProfile[m].nGetY() << "] ";
         //          LogStream << endl;

         double dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

         // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point TODO is this correct? Should it be dParProfStartElev?
         double dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev();

         // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

         double dParProfTotDiff = 0;
         for (int m = 0; m < nParProfLen; m++)
         {
            CGeom2DIPoint PtiTmp = PtiVParProfile[m];
            int
                nX = PtiTmp.nGetX(),
                nY = PtiTmp.nGetY();

            // Safety check
            if (!bIsWithinValidGrid(nX, nY))
            {
               //                LogStream << WARN << "10 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (!m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
            {
               double
                   dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                   dDiff = VdParProfileDeanElev[m] - dTmpElev;

               dParProfTotDiff += dDiff;
            }
         }

         //          // DEBUG STUFF -----------------------------------------------------
         //          LogStream << endl << "\tFor polygon " << nPoly << " doing DOWN-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nUpCoastDeanLen = " << nUpCoastDeanLen << " dUpCoastDeanLen = " << dUpCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
         //
         //          LogStream << "\tExisting profile for deposition = ";
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             CGeom2DIPoint PtiTmp = PtiVParProfile[n];
         //             int
         //                nX = PtiTmp.nGetX(),
         //                nY = PtiTmp.nGetY();
         //
         //             // Safety check
         //             if (! bIsWithinValidGrid(nX, nY))
         //                KeepWithinValidGrid(nX, nY);
         //
         //             LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
         //          }
         //          LogStream << endl;
         //          LogStream << "\tParallel Dean equilibrium profile for deposition = ";
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             LogStream << VdParProfileDeanElev[n] << " ";
         //          }
         //          LogStream << endl;
         //
         //          LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
         //          for (int n = 0; n < nParProfLen; n++)
         //          {
         //             CGeom2DIPoint PtiTmp = PtiVParProfile[n];
         //             int
         //                nX = PtiTmp.nGetX(),
         //                nY = PtiTmp.nGetY();
         //
         //             // Safety check
         //             if (! bIsWithinValidGrid(nX, nY))
         //                KeepWithinValidGrid(nX, nY);
         //
         //             double
         //                dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
         //                dDiff = VdParProfileDeanElev[n] - dTmpElev;
         //
         //             LogStream << dDiff << " ";
         //          }
         //          LogStream << endl;
         //          // END DEBUG STUFF -----------------------------------------------------

         // So will we be able to deposit as much as is needed?
         if (dParProfTotDiff >= dAllSedimentTargetPerProfile)
         {
            //             LogStream << "\tDOWN-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE DEPOSITED" << endl;

            break;
         }
      }

      //      assert(dParProfTotDiff > 0);

      // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the landward end
      double
          dSandDeposited = 0, // Totals for this parallel profile
          dCoarseDeposited = 0,
          dSandToDepositOnProf = dSandTargetPerProfile,
          dCoarseToDepositOnProf = dCoarseTargetPerProfile,
          dSandRatio = dSandToDepositOnProf / (dSandToDepositOnProf + dCoarseToDepositOnProf),
          dCoarseRatio = 1 - dSandRatio;
      //       assert(dRatio >= 0);

      for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
      {
         // Don't bother with tiny amounts
         if (dSandToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
            dSandToDepositOnPoly = 0;
         if (dSandToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
            dSandToDepositOnProf = 0;
         if (dCoarseToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
            dCoarseToDepositOnPoly = 0;
         if (dCoarseToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
            dCoarseToDepositOnProf = 0;

         // Leave the loop if we have deposited enough for this polygon
         if ((dSandToDepositOnPoly <= 0) && (dCoarseToDepositOnPoly <= 0))
         {
            if (dSandToDepositOnPoly < 0)
            {
               LogStream << "!!!!!! DOWN-COAST dSandToDepositOnPoly = " << dSandToDepositOnPoly << endl;
               
               LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
            }

            if (dCoarseToDepositOnPoly < 0)
            {
               LogStream << "!!!!!! DOWN-COAST dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
               
               LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
            }
            
            break;
         }

         // Leave the loop if we have done enough deposition for this profile
         if ((dSandToDepositOnProf <= 0) && (dCoarseToDepositOnProf <= 0))
         {
            if (dSandToDepositOnProf < 0)
            {
               LogStream << "***** DOWN-COAST dSandToDepositOnProf = " << dSandToDepositOnProf << endl;
               
               LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
            }
            
            if (dCoarseToDepositOnProf < 0)
            {
               LogStream << "***** DOWN-COAST dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << endl;
               
               LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
            }
            
            break;
         }

         //          assert(nSeawardFromCoast < PtiVParProfile.size());
         CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
         int
             nX = PtiTmp.nGetX(),
             nY = PtiTmp.nGetY();

         // Safety check
         if (!bIsWithinValidGrid(nX, nY))
         {
            //             LogStream << WARN << "11 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
            PtiTmp.SetX(nX);
            PtiTmp.SetY(nY);
         }

         // Don't do anything to intervention cells
         if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
            continue;

         // Don't do cells twice
         if (!m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
         {
            // Get this cell's current elevation
            double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

            //             LogStream << "\tnPoly = " << nPoly << ", [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

            // Subtract the two elevations
            //             assert(nSeawardFromCoast < VdParProfileDeanElev.size());
            double dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
            if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
            {
               bool
                   bSandDeposited = false,
                   bCoarseDeposited = false;
               double
                   dSandToDeposit = 0,
                   dCoarseToDeposit = 0;

               // The current elevation is below the Dean elevation, so we have can have beach deposition here
               if (dSandToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
               {
                  dSandToDeposit = tMin(dElevDiff * dSandRatio, dSandToDepositOnProf, dSandToDepositOnPoly);

                  int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     bSandDeposited = true;

                     double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisIter[nTopLayer] = true;

                     dSandDeposited += dSandToDeposit;
                     dTotSandDeposited += dSandToDeposit;

                     dSandToDepositOnProf -= dSandToDeposit;
                     dSandToDepositOnPoly -= dSandToDeposit;

                     // Update the cell's beach deposition, and total beach deposition, values
                     m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit);

                     //                      LogStream << m_ulIter << ": nPoly = " << nPoly << ", sand deposition = " << dSandToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                  }
               }

               if (dCoarseToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
               {
                  dCoarseToDeposit = tMin(dElevDiff * dCoarseRatio, dCoarseToDepositOnProf, dCoarseToDepositOnPoly);

                  int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     bCoarseDeposited = true;

                     double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisIter[nTopLayer] = true;

                     dCoarseDeposited += dCoarseToDeposit;
                     dTotCoarseDeposited += dCoarseToDeposit;

                     dCoarseToDepositOnProf -= dCoarseToDeposit;
                     dCoarseToDepositOnPoly -= dCoarseToDeposit;

                     // Update the cell's beach deposition, and total beach deposition, values
                     m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dCoarseToDeposit);

                     //                   LogStream << m_ulIter << ": nPoly = " << nPoly << ", coarse deposition = " << dCoarseToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                  }
               }

               if (bSandDeposited || bCoarseDeposited)
               {
                  // Update the cell's layer elevations
                  m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                  // Update the cell's sea depth
                  m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                  // And set the landform category
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                  // Update this-timestep totals
                  m_ulThisIterNumBeachDepositionCells++;
                  m_dThisIterBeachDepositionSand += dSandToDeposit;
                  m_dThisIterBeachDepositionCoarse += dCoarseToDeposit;
               }
            }
            else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
            {
               // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
               m_ulThisIterNumPotentialBeachErosionCells++;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

               //                   LogStream << m_ulIter << ": nPoly = " << nPoly << " going UP-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Now get the number of the highest layer with non-zero thickness
               int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               // Safety check
               if (nThisLayer == INT_NODATA)
                  return RTN_ERR_NO_TOP_LAYER;

               if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
               {
                  // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                  double
                      dFine = 0,
                      dSand = 0,
                      dCoarse = 0;

                  ErodeBeachSedimentOnCellSupplyLimited(nX, nY, nThisLayer, -dElevDiff, dFine, dSand, dCoarse);

                  // TODO what if dFine is non zero?
                  double dTmpTot = dFine + dSand + dCoarse;
                  if (dTmpTot > 0)
                  {
                     // Update totals for this parallel profile
                     dSandDeposited -= dSand;
                     dCoarseDeposited -= dCoarse;

                     //                      assert(dSandDeposited >= 0);
                     //                      assert(dCoarseDeposited >= 0);

                     // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                     dTotSandDeposited -= dSand;
                     dTotCoarseDeposited -= dCoarse;

                     dSandToDepositOnProf += dSand;
                     dSandToDepositOnPoly += dSand;
                     dCoarseToDepositOnProf += dCoarse;
                     dCoarseToDepositOnPoly += dCoarse;

                     // Update this-timestep totals
                     m_ulThisIterNumActualBeachErosionCells++;
                     m_dThisIterActualBeachErosionFine += dFine;
                     m_dThisIterActualBeachErosionSand += dSand;
                     m_dThisIterActualBeachErosionCoarse += dCoarse;

                     //                      LogStream << m_ulIter << ": in polygon " << nPoly << " have net deposition, actual beach erosion = " << dFine + dSand + dCoarse << " at [" << nX << "][" << nY << "]" << endl;
                  }
               }
            }
         }
      }

      //       LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
   }

   // Have we deposited the full amount?
   double dTotDeposited = dTotSandDeposited + dTotCoarseDeposited;
   if ((!bFPIsEqual(dTotDeposited, dTotTargetToDeposit, TOLERANCE)) && (dTotDeposited < dTotTargetToDeposit))
   {
      // No, so do the same in an UP-COAST (i.e. decreasing coastpoint indices) direction
      // UP-COAST (decreasing coastpoint indices) ==================================================================================

      // Start by finding the seaward end point of the down-coast part-profile, as before we are using only part of each profile, seaward as far as the depth of closure
      //       CGeom2DIPoint PtiDownCoastPartProfileSeawardEnd;
      //       int nIndex1 = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiDownCoastPartProfileSeawardEnd);
      int nIndex1 = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
      if (nIndex1 == INT_NODATA)
      {
         LogStream << m_ulIter << ": " << ERR << "while depositing beach for coast " << nCoast << " polygon " << nPoly << ", could not find the seaward end point of the down-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << endl;
         return RTN_ERR_NO_SEAWARD_END_OF_PROFILE_4;
      }

      // The part-profile length is one greater than nIndex1, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
      int nDownCoastDeanLen = nIndex1 + 1;

      //       assert(bIsWithinValidGrid(&PtiDownCoastPartProfileSeawardEnd));

      // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
      //       CGeom2DPoint
      //          PtDownCoastProfileStart = *pDownCoastProfile->pPtGetPointInProfile(0),
      //          PtDownCoastProfileEnd = PtGridCentroidToExt(&PtiDownCoastPartProfileSeawardEnd);

      //       double dDownCoastDeanLen = dGetDistanceBetween(&PtDownCoastProfileStart, &PtDownCoastProfileEnd);

      int
          nXDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetX(),
          nYDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetY();

      // int nCoastSegLen;

      // Store the coast point numbers for this polygon so that we can shuffle them
      nVCoastPoint.resize(0);
      if (nUpCoastProfileCoastPoint == 0)
      {
         // This is the final up-coast polygon, so also include the up-coast polygon boundary
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

      // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
      Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

      // Recalc the targets for deposition per profile
      dSandTargetPerProfile = dSandToDepositOnPoly / nCoastSegLen;
      dCoarseTargetPerProfile = dCoarseToDepositOnPoly / nCoastSegLen;

      // Now traverse the polygon's existing coastline in a random (but broadly UP-COAST, i.e. decreasing coastpoint indices) sequence, fitting a Dean profile at each coast point
      for (int n = 0; n < nCoastSegLen; n++)
      {
         // Pick a random coast point
         int nCoastPoint = nVCoastPoint[n];
         CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         int
             nCoastX = PtiCoastPoint.nGetX(),
             nCoastY = PtiCoastPoint.nGetY();

         // Calculate the x-y offset between this coast point, and the coast point of the down-coast normal
         int
             nXOffset = nCoastX - nXDownCoastProfileExistingCoastPoint,
             nYOffset = nCoastY - nYDownCoastProfileExistingCoastPoint,
             nSeawardOffset = -1;
         //             nParProfLen;
         vector<CGeom2DIPoint> PtiVParProfile;
         vector<double> VdParProfileDeanElev;

         // OK, loop until we can deposit sufficient unconsolidated sediment
         while (true)
         {
            // Move seaward by one cell
            nSeawardOffset++;

            // And lengthen the parallel profile
            int nParProfLen = nDownCoastDeanLen + nSeawardOffset;

            if (nParProfLen > (pDownCoastProfile->nGetNumCellsInProfile()))
            {
               // We've reached the seaward end of the down-coast profile, need to quit
               // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               //    LogStream << m_ulIter << ": " << WARN << "reached seaward end of down-coast profile during UP-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

               break;
            }

            // Get the x-y coords of a profile starting from this coast point and parallel to the down-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
            PtiVParProfile.resize(0);
            for (int m = 0; m < nParProfLen; m++)
            {
               CGeom2DIPoint PtiProf = *pDownCoastProfile->pPtiGetCellInProfile(m);
               CGeom2DIPoint PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
               PtiVParProfile.push_back(PtiTmp);
            }

            // Get the existing elevation of the seaward end of the parallel profile
            int
                nSeaEndX = PtiVParProfile.back().nGetX(),
                nSeaEndY = PtiVParProfile.back().nGetY();

            // Safety check
            if (!bIsWithinValidGrid(nSeaEndX, nSeaEndY))
            {
               //                LogStream << WARN << "12 @@@@ while doing UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
               PtiVParProfile.back().SetX(nSeaEndX);
               PtiVParProfile.back().SetY(nSeaEndY);
            }

            double dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetSedimentTopElev();

            // Set the start elevation for the Dean profile just a bit above SWL for this timestep, so that it is a Bruun profile
            double dParProfStartElev = m_dThisIterSWL + m_dDeanProfileStartAboveSWL;

            // Calculate the total length of the parallel profile, including any seaward offset
            double dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

            // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
            double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

            // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
            double dParProfA = (dParProfStartElev - dParProfEndElev) / pow(dParProfDeanLen, DEAN_POWER);

            nParProfLen = static_cast<int>(PtiVParProfile.size());
            VdParProfileDeanElev.resize(nParProfLen, 0);

            double dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

            // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point TODO is this correct? Should it be dParProfStartElev?
            double dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev();

            // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
            CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

            double dParProfTotDiff = 0;
            for (int m = 0; m < nParProfLen; m++)
            {
               CGeom2DIPoint PtiTmp = PtiVParProfile[m];
               int
                   nX = PtiTmp.nGetX(),
                   nY = PtiTmp.nGetY();

               // Safety check
               if (!bIsWithinValidGrid(nX, nY))
               {
                  //                   LogStream << WARN << "13 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

                  KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
                  PtiTmp.SetX(nX);
                  PtiTmp.SetY(nY);
               }

               // Don't do anything to intervention cells
               if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
                  continue;

               // Don't do cells twice
               if (!m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
               {
                  double
                      dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                      dDiff = VdParProfileDeanElev[m] - dTmpElev;

                  dParProfTotDiff += dDiff;
               }
            }

            //             // DEBUG STUFF -----------------------------------------------------
            //             LogStream << endl << "\tFor polygon " << nPoly << " doing UP-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nDownCoastDeanLen = " << nDownCoastDeanLen << " dDownCoastDeanLen = " << dDownCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
            //
            //             LogStream << "\tExisting profile for deposition = ";
            //             for (int n = 0; n < nParProfLen; n++)
            //             {
            //                CGeom2DIPoint PtiTmp = PtiVParProfile[n];
            //                int
            //                   nX = PtiTmp.nGetX(),
            //                   nY = PtiTmp.nGetY();
            //
            //                // Safety check
            //                if (! bIsWithinValidGrid(nX, nY))
            //                   KeepWithinValidGrid(nX, nY);
            //
            //                LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
            //             }
            //             LogStream << endl;
            //             LogStream << "\tParallel Dean equilibrium profile for deposition = ";
            //             for (int n = 0; n < nParProfLen; n++)
            //             {
            //                LogStream << VdParProfileDeanElev[n] << " ";
            //             }
            //             LogStream << endl;
            //
            //             LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
            //             for (int n = 0; n < nParProfLen; n++)
            //             {
            //                CGeom2DIPoint PtiTmp = PtiVParProfile[n];
            //                int
            //                   nX = PtiTmp.nGetX(),
            //                   nY = PtiTmp.nGetY();
            //
            //                // Safety check
            //                if (! bIsWithinValidGrid(nX, nY))
            //                   KeepWithinValidGrid(nX, nY);
            //
            //                double
            //                   dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
            //                   dDiff = VdParProfileDeanElev[n] - dTmpElev;
            //
            //                LogStream << dDiff << " ";
            //             }
            //             LogStream << endl;
            //             // END DEBUG STUFF -----------------------------------------------------

            // So will we be able to deposit as much as is needed?
            if (dParProfTotDiff >= dAllSedimentTargetPerProfile)
            {
               //                LogStream << "\tUP-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE DEPOSITED" << endl;

               break;
            }
         }

         //      assert(dParProfTotDiff > 0);

         // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the landward end
         double
             //            dRatio = tMin(1.0, dAllSedimentTargetPerProfile / dParProfTotDiff),
             dSandDeposited = 0, // Totals for this parallel profile
             dCoarseDeposited = 0,
             dSandToDepositOnProf = dSandTargetPerProfile,
             dCoarseToDepositOnProf = dCoarseTargetPerProfile,
             dSandRatio = dSandToDepositOnProf / (dSandToDepositOnProf + dCoarseToDepositOnProf),
             dCoarseRatio = 1 - dSandRatio;

         for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
         {
            // Don't bother with tiny amounts
            if (dSandToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
               dSandToDepositOnPoly = 0;
            if (dSandToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
               dSandToDepositOnProf = 0;
            if (dCoarseToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
               dCoarseToDepositOnPoly = 0;
            if (dCoarseToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
               dCoarseToDepositOnProf = 0;

            // Leave the loop if we have deposited enough for this polygon
            if ((dSandToDepositOnPoly <= 0) && (dCoarseToDepositOnPoly <= 0))
            {
               if (dSandToDepositOnPoly < 0)
               {
                  LogStream << "!!!!!! UP-COAST dSandToDepositOnPoly = " << dSandToDepositOnPoly << endl;
                  
                  LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
               }

               if (dCoarseToDepositOnPoly < 0)
               {
                  LogStream << "!!!!!! UP-COAST dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
                  
                  LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
               }

               break;
            }

            // Leave the loop if we have done enough deposition for this profile
            if ((dSandToDepositOnProf <= 0) && (dCoarseToDepositOnProf <= 0))
            {
               if (dSandToDepositOnProf < 0)
               {
                  LogStream << "***** UP-COAST dSandToDepositOnProf = " << dSandToDepositOnProf << endl;
                  
                  LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
               }
               
               if (dCoarseToDepositOnProf < 0)
               {
                  LogStream << "***** UP-COAST dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << endl;
                  
                  LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
               }
               
               break;
            }

            CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
            int
                nX = PtiTmp.nGetX(),
                nY = PtiTmp.nGetY();

            // Safety check
            if (!bIsWithinValidGrid(nX, nY))
            {
               //                LogStream << WARN << "14 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (!m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisIter())
            {
               // Get this cell's current elevation
               double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

               //          LogStream << "\tnPoly = " << nPoly << " going UP-COAST, [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Subtract the two elevations
               //                assert(nSeawardFromCoast < VdParProfileDeanElev.size());
               double dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
               if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
               {
                  bool
                      bSandDeposited = false,
                      bCoarseDeposited = false;
                  double
                      dSandToDeposit = 0,
                      dCoarseToDeposit = 0;

                  // The current elevation is below the Dean elevation, so we have can have beach deposition here
                  if (dSandToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dSandToDeposit = tMin(dElevDiff * dSandRatio, dSandToDepositOnProf, dSandToDepositOnPoly);

                     int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER;

                     if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        bSandDeposited = true;

                        double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        dSandDeposited += dSandToDeposit;
                        dTotSandDeposited += dSandToDeposit;

                        dSandToDepositOnProf -= dSandToDeposit;
                        dSandToDepositOnPoly -= dSandToDeposit;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit);

                        //                         LogStream << "\tIn nBeachRedistributionOnPolygon() going UP-COAST, nPoly = " << nPoly << " going up-coast, sand deposition = " << dSandToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                  }

                  if (dCoarseToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dCoarseToDeposit = tMin(dElevDiff * dCoarseRatio, dCoarseToDepositOnProf, dCoarseToDepositOnPoly);

                     int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER;

                     if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        bCoarseDeposited = true;

                        double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisIter[nTopLayer] = true;

                        dCoarseDeposited += dCoarseToDeposit;
                        dTotCoarseDeposited += dCoarseToDeposit;

                        dCoarseToDepositOnProf -= dCoarseToDeposit;
                        dCoarseToDepositOnPoly -= dCoarseToDeposit;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dCoarseToDeposit);

                        //                         LogStream << m_ulIter << ": nPoly = " << nPoly << " going UP-COAST, coarse deposition = " << dCoarseToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                  }

                  if (bSandDeposited || bCoarseDeposited)
                  {
                     // Now update the cell's layer elevations
                     m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                     // Update the cell's sea depth
                     m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                     // And set the landform category
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                     // Update this-timestep totals
                     m_ulThisIterNumBeachDepositionCells++;
                     m_dThisIterBeachDepositionSand += dSandToDeposit;
                     m_dThisIterBeachDepositionCoarse += dCoarseToDeposit;
                  }
               }
               else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
               {
                  // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
                  m_ulThisIterNumPotentialBeachErosionCells++;

                  m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

                  //                   LogStream << m_ulIter << ": nPoly = " << nPoly << " going UP-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

                  // Now get the number of the highest layer with non-zero thickness
                  int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

                  // Safety check
                  if (nThisLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
                  {
                     // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                     double
                         dFine = 0,
                         dSand = 0,
                         dCoarse = 0;

                     ErodeBeachSedimentOnCellSupplyLimited(nX, nY, nThisLayer, -dElevDiff, dFine, dSand, dCoarse);

                     // TODO what if dFine is non zero?
                     double dTmpTot = dFine + dSand + dCoarse;
                     if (dTmpTot > 0)
                     {
                        // Update totals for this parallel profile
                        dSandDeposited -= dSand;
                        dCoarseDeposited -= dCoarse;

                        //                      assert(dSandDeposited >= 0);
                        //                      assert(dCoarseDeposited >= 0);

                        // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                        dTotSandDeposited -= dSand;
                        dTotCoarseDeposited -= dCoarse;

                        dSandToDepositOnProf += dSand;
                        dSandToDepositOnPoly += dSand;
                        dCoarseToDepositOnProf += dCoarse;
                        dCoarseToDepositOnPoly += dCoarse;

                        // Update this-timestep totals
                        m_ulThisIterNumActualBeachErosionCells++;
                        m_dThisIterActualBeachErosionFine += dFine;
                        m_dThisIterActualBeachErosionSand += dSand;
                        m_dThisIterActualBeachErosionCoarse += dCoarse;

                        //                         LogStream << m_ulIter << ": in polygon " << nPoly << " have net deposition, actual beach erosion = " << dFine + dSand + dCoarse << " at [" << nX << "][" << nY << "]" << endl;
                     }
                  }
               }
            }
         }

         //          LogStream << m_ulIter << ": nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
      }
   }

   // How much have we been able to deposit on this polygon?
   double dTmpTot = dTotSandDeposited + dTotCoarseDeposited;
   string strMsg;
   if (bFPIsEqual(dTmpTot, dTotTargetToDeposit, SEDIMENT_ELEV_TOLERANCE))
      strMsg = "deposition OK";
   else
   {
      if (dTmpTot > dTotTargetToDeposit)
         strMsg = "TOO MUCH deposition";
      else
         strMsg = "NOT ENOUGH deposition";

      if (dTmpTot != 0)
      {
         // m_dThisIterDepositionSandDiff and m_dThisIterDepositionCoarseDiff are both +ve if there is not enough deposition
         m_dThisIterDepositionSandDiff += (dTotTargetToDeposit - dTmpTot) * dTotSandDeposited / dTmpTot;
         m_dThisIterDepositionCoarseDiff += (dTotTargetToDeposit - dTmpTot) * dTotCoarseDeposited / dTmpTot;
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": polygon " << nPoly << " " << strMsg << " deposited = " << (dTotSandDeposited + dTotCoarseDeposited) * m_dCellArea << " target = " << dTotTargetToDeposit * m_dCellArea << " (sand = " << dTotSandDeposited * m_dCellArea << " sand remaining = " << dSandToDepositOnPoly * m_dCellArea << " coarse = " << dTotCoarseDeposited * m_dCellArea << " coarse remaining = " << dCoarseToDepositOnPoly * m_dCellArea << "), all m^3" << endl;

   return RTN_OK;
}
