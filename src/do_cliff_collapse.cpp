/*!
 *
 * \file do_cliff_collapse.cpp
 * \brief Collapses cliffs
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

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include "cme.h"
#include "simulation.h"
#include "cliff.h"

/*===============================================================================================================================

 Update accumulated wave energy in coastal landform objects

===============================================================================================================================*/
int CSimulation::nDoAllWaveEnergyToCoastLandforms(void)
{
   int nRet;

   // First go along each coastline and update the total wave energy which it has experienced
   for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
   {
      for (int j = 0; j < m_VCoast[i].nGetCoastlineSize(); j++)
      {
         CACoastLandform *pCoastLandform = m_VCoast[i].pGetCoastLandform(j);

         // Update accumulated wave energy for the coastal landform object
         double
             dWaveHeightAtCoast = m_VCoast[i].dGetCoastWaveHeight(j),
             dDeepWaterWavePeriod = m_VCoast[i].dGetCoastDeepWaterWavePeriod(j),
             dWaveErosiveForce = pow(dWaveHeightAtCoast, WALKDEN_HALL_PARAM_1) * pow(dDeepWaterWavePeriod, WALKDEN_HALL_PARAM_2),
             dWaveEnergy = dWaveErosiveForce * m_dTimeStep / 3600 / 1000;

         //          assert(isfinite(dWaveEnergy));
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         // Now simulate how the coastal landform responds to this wave energy
         int nCategory = pCoastLandform->nGetLandFormCategory();
         if (nCategory == LF_CAT_CLIFF)
         {
            // This is a cliff
            CRWCliff *pCliff = reinterpret_cast<CRWCliff *>(pCoastLandform);

            // Calculate this-timestep cliff notch erosion (is a length in external CRS units)
            double dNotchExtension = dWaveEnergy / m_dCliffErosionResistance;

            // Constrain this-timestep notch extension if it is more than the length of one cell side (in external CRS units), since the most we can remove in a single timestep is one coastal cell
            dNotchExtension = tMin(m_dCellSide, dNotchExtension);

            // Extend the cliff object's erosional notch as a result of wave energy during this timestep. Note that extension may be constrained, since this-timestep extension cannot exceed the depth of sediment remaining on the cell
            dNotchExtension = pCliff->dErodeNotch(dNotchExtension); // Constrain

            // OK, is the notch now extended enough to cause collapse (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchOverhangAtCollapse))
            {
               // It is ready to collapse
               double
                   dFineCollapse = 0,
                   dSandCollapse = 0,
                   dCoarseCollapse = 0;

               // So do the cliff collapse
               nRet = nDoCliffCollapse(pCliff, dNotchExtension, dFineCollapse, dSandCollapse, dCoarseCollapse);
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     LogStream << m_ulIter << WARN << " problem with cliff collapse, continuing however" << endl;
               }

               // And put fine sediment into suspension, and deposit sand and/or coarse sediment as unconsolidated sediment
               nRet = nDoCliffCollapseDeposition(pCliff, dFineCollapse, dSandCollapse, dCoarseCollapse);
               if (nRet != RTN_OK)
                  return nRet;
            }
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << endl
                << m_ulIter << ": cliff collapse (m^3) = " << (m_dThisIterCliffCollapseErosionFine + m_dThisIterCliffCollapseErosionSand + m_dThisIterCliffCollapseErosionCoarse) * m_dCellArea << " (fine = " << m_dThisIterCliffCollapseErosionFine * m_dCellArea << ", sand = " << m_dThisIterCliffCollapseErosionSand * m_dCellArea << ", coarse = " << m_dThisIterCliffCollapseErosionCoarse * m_dCellArea << "), talus deposition (m^3) = " << (m_dThisIterCliffDepositionSand + m_dThisIterCliffDepositionCoarse) * m_dCellArea << " (sand = " << m_dThisIterCliffDepositionSand * m_dCellArea << ", coarse = " << m_dThisIterCliffDepositionSand * m_dCellArea << ")" << endl;

   return RTN_OK;
}

/*===============================================================================================================================

 Simulates cliff collapse on a single cliff object: it updates both the cliff object and the cell 'under' the cliff object

===============================================================================================================================*/
int CSimulation::nDoCliffCollapse(CRWCliff *pCliff, double const dNotchDeepen, double &dFineCollapse, double &dSandCollapse, double &dCoarseCollapse)
{
   // Get the cliff cell's grid coords
   int
       nX = pCliff->pPtiGetCellMarkedAsLF()->nGetX(),
       nY = pCliff->pPtiGetCellMarkedAsLF()->nGetY();

   // Then get the elevation of the base of the notch from the cliff object
   double dNotchElev = pCliff->dGetNotchBaseElev() - m_dNotchBaseBelowSWL;

   // Get the index of the layer containing the notch (layer 0 being just above basement)
   int nNotchLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dNotchElev);
   if (nNotchLayer == ELEV_ABOVE_SEDIMENT_TOP)
   {
      if (m_nLogFileDetail >= LOG_FILE_MOST_DETAIL)
         LogStream << endl
                   << m_ulIter << ": " << ERR << " cell [" << nX << "][" << nY << "] has dNotchElev (" << dNotchElev << ") above sediment top elevation (" << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << ")" << endl;

      return RTN_ERR_CLIFFNOTCH;
   }
   //    else
   //       LogStream << endl << m_ulIter << ": for cell [" << nX << "][" << nY << "] dNotchElev = " << dNotchElev << " sediment top elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << endl;

   // Flag the coastline cliff object as having collapsed
   pCliff->SetCliffCollapse(true);

   int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

   // Safety check
   if (nTopLayer == INT_NODATA)
      return RTN_ERR_NO_TOP_LAYER;

   double dRemaining = pCliff->dGetRemaining();
   if (dRemaining <= 0)
   {
      // No cliff sediment left on this cliff object, so the cell which it occupies will no longer be a cliff in the next timestep
      pCliff->SetAllSedimentGone();

      // Set the base of the collapse (see above)
      pCliff->SetNotchBaseElev(dNotchElev);

      // Set flags to say that the top layer has changed
      m_bConsChangedThisIter[nTopLayer] = true;
      m_bUnconsChangedThisIter[nTopLayer] = true;

      //      int nX = pCliff->pPtiGetCellMarkedAsLF()->nGetX();
      //      int nY = pCliff->pPtiGetCellMarkedAsLF()->nGetY();
      //      LogStream << m_ulIter << ": all sediment removed from cliff object after cliff collapse on [" << nX << "][" << nY << "], dNotchElev = " << dNotchElev << endl;
   }

   // Now calculate the vertical depth of sediment lost in this cliff collapse
   double dAboveNotch = m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() - dNotchElev;

   // In CoastalME, all depth equivalents are assumed to be a depth upon the whole of a cell i.e. upon the area of a whole cell. The vertical depth of sediment lost in each cliff collapse is a depth upon only part of a cell, i.e. upon a fraction of a cell's area. To keep the depth of cliff collapse consistent with all other depth equivalents, weight it by the fraction of the cell's area which is being removed
   double
       dNotchAreaFrac = dNotchDeepen / m_dCellSide,
       dCollapseDepth = dAboveNotch * dNotchAreaFrac;

   // Update the cell's totals for cliff collapse
   m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapse(dCollapseDepth);

   double
       dAvailable = 0,
       dLost;

   //   LogStream << m_ulIter << ": cell [" << nX << "][" << nY << "] before removing sediment, dGetVolEquivSedTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() << ", dGetSedimentTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << endl;

   // Now update the cell's sediment. If there are sediment layers above the notched layer, we must remove sediment from the whole depth of each layer. Again, weight the depth lost by the fraction of the cell's area which is being removed
   for (int n = nTopLayer; n > nNotchLayer; n--)
   {
      // Add in this layer's sediment, both consolidated and unconsolidated, and adjust what is left
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dFineCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchFineLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dSandCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchSandLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dCoarseCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchCoarseLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dFineCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchFineLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dSandCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchSandLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dCoarseCollapse += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchCoarseLost(dLost);
      }
   }

   // For the layer which contains the notch, remove only part of the sediment depth
   double
       dNotchLayerTop = m_pRasterGrid->m_Cell[nX][nY].dCalcLayerElev(nNotchLayer),
       dNotchLayerThickness = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->dGetTotalThickness(),
       dNotchLayerVertFracRemoved = (dNotchLayerTop - dNotchElev) / dNotchLayerThickness;

   // Now calculate the fraction of the volume which is removed
   double dNotchLayerFracRemoved = dNotchLayerVertFracRemoved * dNotchAreaFrac;

   // Sort out the notched layer's sediment, both consolidated and unconsolidated, for this cell
   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchFineLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchSandLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchCoarseLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchFineLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchSandLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchCoarseLost(dLost);
   }

   // Update the cell's layer elevations
   m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

   // And update the cell's sea depth
   m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

   // The notch has gone
   pCliff->SetNotchOverhang(0);

   //   LogStream << m_ulIter << ": cell [" << nX << "][" << nY << "] after removing sediment, dGetVolEquivSedTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() << ", dGetSedimentTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << endl << endl;

   // And update the this-timestep totals and the grand totals for collapse
   m_nNThisIterCliffCollapse++;
   m_nNTotCliffCollapse++;

   m_dThisIterCliffCollapseErosionFine += dFineCollapse;
   m_dThisIterCliffCollapseErosionSand += dSandCollapse;
   m_dThisIterCliffCollapseErosionCoarse += dCoarseCollapse;

   return RTN_OK;
}

/*===============================================================================================================================

 Puts the fine sediment from a cliff collapse into suspension, and redistributes the sand-sized and coarse-sized sediment from a cliff collapse onto the foreshore, as unconsolidated talus

 The talus is added to the existing beach volume (i.e. to the unconsolidated sediment). The shoreline is iteratively advanced seaward until all this volume is accommodated under a Dean equilibrium profile. This equilibrium beach profile is h(y) = A * y^(2/3) where h(y) is the water depth at a distance y from the shoreline and A is a sediment-dependent scale parameter

===============================================================================================================================*/
int CSimulation::nDoCliffCollapseDeposition(CRWCliff *pCliff, double const dFineCollapse, double const dSandCollapse, double const dCoarseCollapse)
{
   // Fine sediment goes into suspension
   if (dFineCollapse > SEDIMENT_ELEV_TOLERANCE)
   {
      m_dThisIterFineSedimentToSuspension += dFineCollapse;
   }

   // Do we have any sand- or coarse-sized sediment to deposit?
   double dTotFromCollapse = dSandCollapse + dCoarseCollapse;
   if (dTotFromCollapse < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   // OK, we have some sand- and/or coarse-sized sediment to deposit
   int
       nCoast = pCliff->nGetCoast(),
       nStartPoint = pCliff->nGetPointOnCoast(),
       nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   double
       dTotSandToDeposit = dSandCollapse,
       dTotCoarseToDeposit = dCoarseCollapse,
       dSandProp = dSandCollapse / dTotFromCollapse,
       dCoarseProp = 1 - dSandProp;

   //    LogStream << "=====================================================================================================" << endl;
   //    LogStream << m_ulIter << ": coast = " << nCoast << ", point = " << nStartPoint << endl;

   // Calculate the proportion per planview collapse profile
   vector<int>
       nVWidthDistSigned(m_nCliffCollapseTalusPlanviewWidth),
       nVProfileLength(m_nCliffCollapseTalusPlanviewWidth);
   vector<double>
       dVToDepositPerProfile(m_nCliffCollapseTalusPlanviewWidth);

   // TODO IMPROVE LATER, at present all planview profiles are the same length
   int nSigned = -(m_nCliffCollapseTalusPlanviewWidth - 1) / 2;
   for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
   {
      nVWidthDistSigned[n] = nSigned++;
      nVProfileLength[n] = nRound(m_dCliffDepositionPlanviewLength);
      dVToDepositPerProfile[n] = (dTotSandToDeposit + dTotCoarseToDeposit) / m_nCliffCollapseTalusPlanviewWidth;
   }

   //    LogStream << "Width offsets = ";
   //    for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
   //    {
   //       LogStream << nVWidthDistSigned[n] << " ";
   //    }
   //    LogStream << endl << "Profile lengths = ";
   //    for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
   //    {
   //       LogStream << nVProfileLength[n] << " ";
   //    }
   //    LogStream << endl << "Deposition per profile = ";
   //    for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
   //    {
   //       LogStream << dVToDepositPerProfile[n] << " ";
   //    }
   //    LogStream << endl;

   //   LogStream << "dSandCollapse = " << dSandCollapse << " dCoarseCollapse = " << dCoarseCollapse << " m_nCliffCollapseTalusPlanviewWidth = " << m_nCliffCollapseTalusPlanviewWidth << endl;

   //    int nX = pCliff->pPtiGetCellMarkedAsLF()->nGetX();
   //    int nY = pCliff->pPtiGetCellMarkedAsLF()->nGetY();
   //   LogStream << "Cliff object is at cell[" << nX << "][" << nY << "] which is " << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << endl;

   for (int nAcross = 0; nAcross < m_nCliffCollapseTalusPlanviewWidth; nAcross++)
   {
      int nWidthDistSigned = nVWidthDistSigned[nAcross];

      // Get the start point of the deposition collapse profile
      int nThisPoint = nStartPoint + nWidthDistSigned;

      // Is this start point valid?
      if ((nThisPoint < 0) || (nThisPoint > (nCoastSize - 1)))
      {
         //          LogStream << endl << m_ulIter << ": ABANDONING PROFILE with nWidthDistSigned = " << nWidthDistSigned << endl;
         //          LogStream << "START point " << nThisPoint << " of profile would have been outside the grid, so " << dVToDepositPerProfile[nAcross] << " exported from grid" << endl;
         //          LogStream << "dTotSandToDeposit WAS = " << dTotSandToDeposit << " dTotCoarseToDeposit WAS = " << dTotCoarseToDeposit << endl;

         // The start point of the profile would have been outside the grid, so just add this profile's sediment to the volume exported from the grid this timestep
         m_dThisIterSandSedLostCliffCollapse += (dVToDepositPerProfile[nAcross] * dSandProp);
         m_dThisIterCoarseSedLostCliffCollapse += (dVToDepositPerProfile[nAcross] * dCoarseProp);

         // Remove this volume from the total still to be deposited
         dTotSandToDeposit -= (dVToDepositPerProfile[nAcross] * dSandProp);
         dTotCoarseToDeposit -= (dVToDepositPerProfile[nAcross] * dCoarseProp);
         //          LogStream << "dTotSandToDeposit NOW = " << dTotSandToDeposit << " dTotCoarseToDeposit NOW = " << dTotCoarseToDeposit << endl;

         continue;
      }

      CGeom2DPoint
          PtStart,
          PtEnd;

      // Make the start of the deposition profile the cliff cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
      PtStart.SetX(dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetX()));
      PtStart.SetY(dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetY()));

      // The seaward offset, in cells
      int nSeawardOffset = -1;
      do
      {
         //         LogStream << endl;

         nSeawardOffset++;

         //          if (nSeawardOffset > 20)
         //          {
         //             // Arbitrary safety check, if we can't store sufficient sediment with an offset of this size then move to the next point along the coast
         //             LogStream << "*** PROFILE TOO LONG WITH nSeawardOffset = " << nSeawardOffset << " MOVING TO NEXT POINT ALONG COAST" << endl;
         //             break;
         //          }

         // Now construct a deposition collapse profile from the start point, it is one longer than the specified length because it includes the cliff point in the profile. Calculate its length in external CRS units, the way it is done here is approximate but probably OK
         double dThisProfileLength = (nVProfileLength[nAcross] + nSeawardOffset + 1) * m_dCellSide;

         // Get the end point of this coastline-normal line
         CGeom2DIPoint PtiEnd; // In grid CRS
         int nRtn = nGetCoastNormalEndPoint(nCoast, nThisPoint, nCoastSize, &PtStart, dThisProfileLength, &PtEnd, &PtiEnd);
         if (nRtn != RTN_OK)
         {
            // Could not find an end point so forget this profile
            //            LogStream << endl << m_ulIter << ": ABANDONING PROFILE with nWidthDistSigned = " << nWidthDistSigned << endl;

            if (nRtn == RTN_ERR_PROFILE_ENDPOINT_IS_OFFGRID)
            {
               //                LogStream << "dTotSandToDeposit WAS = " << dTotSandToDeposit << " dTotCoarseToDeposit WAS = " << dTotCoarseToDeposit << endl;
               //                LogStream << "END point of profile would have been outside the grid, so " << dVToDepositPerProfile[nAcross] << " exported from grid" << endl;

               // The end point of the profile would have been outside the grid, so just add this profile's sediment to the volume exported from the grid this timestep
               m_dThisIterSandSedLostCliffCollapse += (dVToDepositPerProfile[nAcross] * dSandProp);
               m_dThisIterCoarseSedLostCliffCollapse += (dVToDepositPerProfile[nAcross] * dCoarseProp);

               // Remove this volume from the total still to be deposited
               dTotSandToDeposit -= (dVToDepositPerProfile[nAcross] * dSandProp);
               dTotCoarseToDeposit -= (dVToDepositPerProfile[nAcross] * dCoarseProp);
               //                LogStream << "dTotSandToDeposit NOW = " << dTotSandToDeposit << " dTotCoarseToDeposit NOW = " << dTotCoarseToDeposit << endl;
            }

            if (nRtn == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
            {
               // The profile has a different problem, so (if possible) move one point along the coast and try again. Must add this profile's sediment to the amount remaining per profile however
               //                for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
               //                {
               //                   LogStream << dVToDepositPerProfile[n] << " ";
               //                }
               //                LogStream << endl;
               //                LogStream << "Deposition per profile WAS = ";

               int nWidthRemaining = m_nCliffCollapseTalusPlanviewWidth - nAcross - 1;
               for (int n = nAcross + 1; n < m_nCliffCollapseTalusPlanviewWidth; n++)
                  dVToDepositPerProfile[n] = (dTotSandToDeposit + dTotCoarseToDeposit) / nWidthRemaining;

               //                LogStream << "Deposition per profile NOW = ";
               //                for (int n = 0; n < m_nCliffCollapseTalusPlanviewWidth; n++)
               //                {
               //                   LogStream << dVToDepositPerProfile[n] << " ";
               //                }
               //                LogStream << endl;
            }

            break;
         }

         //         LogStream << m_ulIter << ": nWidthDistSigned = " << nWidthDistSigned << " cliff collapse profile from " << PtStart.dGetX() << ", " << PtStart.dGetY() << " to " << PtEnd.dGetX() << ", " << PtEnd.dGetY() << " with length (inc. cliff point) = " << dThisProfileLength << endl;

         vector<CGeom2DPoint> VTmpProfile;
         VTmpProfile.push_back(PtStart);
         VTmpProfile.push_back(PtEnd);
         vector<CGeom2DIPoint> VCellsUnderProfile;

         // Now get the raster cells under this profile
         if (nRasterizeCliffCollapseProfile(&VTmpProfile, &VCellsUnderProfile) != RTN_OK)
         {
            cout << m_ulIter << ": error when rasterizing cells during cliff collapse" << endl;
            return RTN_ERR_LINETOGRID;
         }

         int nRasterProfileLength = static_cast<int>(VCellsUnderProfile.size());
         vector<double> dVProfileNow(nRasterProfileLength, 0);
         vector<bool> bVProfileValid(nRasterProfileLength, true);
         //         LogStream << "RASTER PROFILE LENGTH = " << nRasterProfileLength << endl;
         //         if (nRasterProfileLength != dThisProfileLength)
         //            LogStream << "*************************" << endl;

         for (int n = 0; n < nRasterProfileLength; n++)
         {
            int
                nX = VCellsUnderProfile[n].nGetX(),
                nY = VCellsUnderProfile[n].nGetY();

            dVProfileNow[n] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

            // Don't allow cliff collapse onto intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               bVProfileValid[n] = false;
         }

         // Calculate the elevation of the talus top
         double
             dCliffTopElev = dVProfileNow[0],
             dCliffBaseElev = dVProfileNow[1],
             dCliffHeight = dCliffTopElev - dCliffBaseElev,
             dTalusTopElev = dCliffBaseElev + (dCliffHeight * m_dCliffDepositionHeightFrac);

         //         LogStream << "Elevations: cliff top = " << dCliffTopElev << " cliff base = " << dCliffBaseElev << " talus top = " << dTalusTopElev << endl;

         //          if (dCliffTopElev < dCliffBaseElev)
         //             LogStream << "*** ERROR, cliff top is lower than cliff base" << endl;

         // The talus slope length in external CRS units, this is approximate but probably OK
         double dTalusSlopeLength = dThisProfileLength - ((nSeawardOffset - 1) * m_dCellSide);

         // If user has not supplied a value for m_dCliffDepositionA, then solve for dA so that the elevations at end of the existing profile, and at the end of the Dean equilibrium profile, are the same
         double dA = 0;
         if (m_dCliffDepositionA != 0)
            dA = m_dCliffDepositionA;
         else
            dA = (dTalusTopElev - dVProfileNow[nRasterProfileLength - 1]) / pow(dTalusSlopeLength, DEAN_POWER);

         double dInc = dTalusSlopeLength / (nRasterProfileLength - nSeawardOffset - 2);
         vector<double> dVDeanProfile(nRasterProfileLength);

         // Calculate the Dean equilibrium profile of the talus h(y) = A * y^(2/3) where h(y) is the distance below the talus-top elevation (the highest point in the Dean profile) at a distance y from the cliff (the landward start of the profile)
         CalcDeanProfile(&dVDeanProfile, dInc, dTalusTopElev, dA, true, nSeawardOffset, dCliffTopElev);

         // Get the total difference in elevation between the two profiles (present profile - Dean profile)
         double dTotElevDiff = dSubtractProfiles(&dVProfileNow, &dVDeanProfile, &bVProfileValid);

         //          // DEBUG STUFF -----------------------------------------------------
         //          LogStream << endl;
         //          LogStream << "dTalusSlopeLength = " << dTalusSlopeLength << " dA = " << dA << endl;
         //          LogStream << "dDistFromTalusStart - dInc = " << dDistFromTalusStart - dInc << " dThisProfileLength - nSeawardOffset - 2 = " << dThisProfileLength - nSeawardOffset - 2 << endl;
         //          LogStream << "Profile now (inc. cliff cell) = ";
         //          for (int n = 0; n < nRasterProfileLength; n++)
         //          {
         //             int
         //                nX = VCellsUnderProfile[n].nGetX(),
         //                nY = VCellsUnderProfile[n].nGetY();
         //             dVProfileNow[n] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();
         //             LogStream << dVProfileNow[n] << " ";
         //          }
         //          LogStream << endl;
         //          LogStream << "Dean equilibrium profile (inc. cliff cell) = ";
         //          for (int n = 0; n < nRasterProfileLength; n++)
         //          {
         //             LogStream << dVDeanProfile[n] << " ";
         //          }
         //          LogStream << endl;
         //          LogStream << "Difference (inc. cliff cell) = ";
         //          for (int n = 0; n < nRasterProfileLength; n++)
         //          {
         //             LogStream << dVDeanProfile[n] - dVProfileNow[n] << " ";
         //          }
         //          LogStream << endl;
         //          // DEBUG STUFF -----------------------------------------------------

         // For this planview profile, does the Dean equilibrium profile allow us to deposit all the talus sediment which we need to get rid of?
         if (dTotElevDiff < dVToDepositPerProfile[nAcross])
            // No it doesn't, so try again with a larger seaward offset
            break;

         // Yes it does
         //          LogStream << m_ulIter << ": cliff collapse at [" << VCellsUnderProfile[0].nGetX() << "][" << VCellsUnderProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(VCellsUnderProfile[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(VCellsUnderProfile[0].nGetY()) << "} offset SUFFICIENT with nSeawardOffset = " << nSeawardOffset << endl;
         //          LogStream << m_ulIter << ": dTotElevDiff = " << dTotElevDiff << " dVToDepositPerProfile[nAcross] = " << dVToDepositPerProfile[nAcross] << endl;

         double dPropToDeposit = dVToDepositPerProfile[nAcross] / dTotElevDiff;
         //         LogStream << "dPropToDeposit = " << dPropToDeposit << endl;

         //          double
         //             dDepositedCheck = 0,
         //             dRemovedCheck = 0;

         // Adjust all cells in this profile
         for (int n = 0; n < nRasterProfileLength; n++)
         {
            int
                nX = VCellsUnderProfile[n].nGetX(),
                nY = VCellsUnderProfile[n].nGetY();

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

            // Safety check
            if (nTopLayer == INT_NODATA)
               return RTN_ERR_NO_TOP_LAYER;

            if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
            {
               // TODO improve this
               cerr << "All layers have zero thickness" << endl;
               return RTN_ERR_CLIFFDEPOSIT;
            }

            // Only do deposition on this cell if its elevation is below the Dean elevation, and the cell is either a sea cell or a drift cell
            if ((dVDeanProfile[n] > dVProfileNow[n]) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
            {
               //               LogStream << "DEPOSIT ";
               // At this point along the profile, the equilibrium profile is higher than the present profile. So we can deposit some sediment here
               double dPotentialSandToDeposit = 0;
               if (dTotSandToDeposit > 0)
               {
                  dPotentialSandToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dSandProp * dPropToDeposit;
                  dPotentialSandToDeposit = tMin(dPotentialSandToDeposit, dTotSandToDeposit);

                  double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dPotentialSandToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotSandToDeposit -= dPotentialSandToDeposit;
                  //                  dDepositedCheck += dPotentialSandToDeposit;
               }

               double dPotentialCoarseToDeposit = 0;
               if (dTotCoarseToDeposit > 0)
               {
                  dPotentialCoarseToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dCoarseProp * dPropToDeposit;
                  dPotentialCoarseToDeposit = tMin(dPotentialCoarseToDeposit, dTotCoarseToDeposit);

                  double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dPotentialCoarseToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotCoarseToDeposit -= dPotentialCoarseToDeposit;
                  //                  dDepositedCheck += dPotentialCoarseToDeposit;
               }

               // Now update the cell's layer elevations
               m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

               // Update the cell's sea depth
               m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

               // Update the cell's collapse deposition, and total collapse deposition, values
               m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapseDeposition(dPotentialSandToDeposit + dPotentialCoarseToDeposit);

               // And set the landform category
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_TALUS);
            }

            else if (dVDeanProfile[n] < dVProfileNow[n])
            {
               // The Dean equilibrium profile is lower than the present profile, so we must remove some some sediment from here
               //               LogStream << "REMOVE ";
               double dThisLowering = dVProfileNow[n] - dVDeanProfile[n];

               // Find out how much sediment we have available on this cell
               double dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetFine();
               double dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
               double dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();

               // Now partition the total lowering for this cell between the three size fractions: do this by relative erodibility
               int nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0);
               int nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0);
               int nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

               double dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);
               //              double dTotActualErosion = 0;

               if (nFineWeight)
               {
                  // Erode some fine-sized sediment
                  double dFineLowering = (m_dFineErodibilityNormalized * dThisLowering) / dTotErodibility;

                  // Make sure we don't get -ve amounts left on the cell
                  double dFine = tMin(dExistingAvailableFine, dFineLowering);
                  double dRemaining = dExistingAvailableFine - dFine;

                  //                 dTotActualErosion += dFine;
                  //                 dRemovedCheck += dFine;

                  // Set the value for this layer
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetFine(dRemaining);

                  // And set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  // And increment the per-timestep talus erosion total
                  m_dThisIterCliffErosionFine += dFine;
               }

               if (nSandWeight)
               {
                  // Erode some sand-sized sediment
                  double dSandLowering = (m_dSandErodibilityNormalized * dThisLowering) / dTotErodibility;

                  // Make sure we don't get -ve amounts left on the source cell
                  double dSand = tMin(dExistingAvailableSand, dSandLowering);
                  double dRemaining = dExistingAvailableSand - dSand;

                  //                 dTotActualErosion += dSand;
                  //                 dRemovedCheck += dSand;

                  // Set the value for this layer
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dRemaining);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  // And increment the per-timestep talus erosion total
                  m_dThisIterCliffTalusSandErosion += dSand;
               }

               if (nCoarseWeight)
               {
                  // Erode some coarse-sized sediment
                  double dCoarseLowering = (m_dCoarseErodibilityNormalized * dThisLowering) / dTotErodibility;

                  // Make sure we don't get -ve amounts left on the source cell
                  double dCoarse = tMin(dExistingAvailableCoarse, dCoarseLowering);
                  double dRemaining = dExistingAvailableCoarse - dCoarse;

                  //                 dTotActualErosion += dCoarse;
                  //                 dRemovedCheck += dCoarse;

                  // Set the value for this layer
                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dRemaining);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  // And increment the per-timestep talus erosion total
                  m_dThisIterCliffTalusCoarseErosion += dCoarse;
               }

               // Set the actual erosion value for this cell
               //               m_pRasterGrid->m_Cell[nX][nY].SetActualPlatformErosion(dTotActualErosion);

               // Recalculate the elevation of every layer
               m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

               // And update the cell's sea depth
               m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

               //                // Update per-timestep totals
               //                if (dTotActualErosion > 0)
               //                {
               //                   m_ulThisIterNumActualPlatformErosionCells++;
               //                   m_dThisIterActualPlatformErosion += dTotActualErosion;
               //                }
            }
         }
         //        LogStream << endl;
         //        LogStream << "Profile done, dDepositedCheck = " << dDepositedCheck << " dRemovedCheck = " << dRemovedCheck << endl;

         break;
      } while (true);

      //     LogStream << "REMAINING dTotSandToDeposit = " << dTotSandToDeposit << " dTotCoarseToDeposit = " << dTotCoarseToDeposit << endl;
   }

   // Safety check
   //   if (! bFPIsEqual(dTotSandToDeposit, 0, TOLERANCE))
   //   {
   //      LogStream << ERR << m_ulIter << ": dTotSandToDeposit = " << dTotSandToDeposit << " SET TO ZERO" << endl;
   //      dTotSandToDeposit = 0;
   //   }

   // Ditto
   //   if (! bFPIsEqual(dTotCoarseToDeposit, 0, TOLERANCE))
   //   {
   //      LogStream << ERR << m_ulIter << ": dTotCoarseToDeposit = " << dTotCoarseToDeposit << " SET TO ZERO" << endl;
   //      dTotCoarseToDeposit = 0;
   //   }

   // Increment this-timestep totals for cliff collapse deposition
   m_dThisIterCliffDepositionSand += dSandCollapse;
   m_dThisIterCliffDepositionCoarse += dCoarseCollapse;

   return RTN_OK;
}

/*==============================================================================================================================

 Given the start and end points of a cliff-collapse normal profile, returns an output vector of cells which are 'under' the vector line

===============================================================================================================================*/
int CSimulation::nRasterizeCliffCollapseProfile(vector<CGeom2DPoint> const *pVPointsIn, vector<CGeom2DIPoint> *pVIPointsOut) const
{
   pVIPointsOut->clear();

   // The start point of the normal is the centroid of a coastline cell. Convert from the external CRS to grid CRS
   double
       dXStart = dExtCRSXToGridX(pVPointsIn->at(0).dGetX()),
       dYStart = dExtCRSYToGridY(pVPointsIn->at(0).dGetY());

   // The end point of the normal, again convert from the external CRS to grid CRS. Note too that it could be off the grid
   double
       dXEnd = dExtCRSXToGridX(pVPointsIn->at(1).dGetX()),
       dYEnd = dExtCRSYToGridY(pVPointsIn->at(1).dGetY());

   // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
   double
       dXInc = dXEnd - dXStart,
       dYInc = dYEnd - dYStart,
       dLength = tMax(tAbs(dXInc), tAbs(dYInc));

   dXInc /= dLength;
   dYInc /= dLength;

   double
       dX = dXStart,
       dY = dYStart;

   // Process each interpolated point
   int nLength = nRound(dLength);
   for (int m = 0; m <= nLength; m++)
   {
      int
          nX = static_cast<int>(dX),
          nY = static_cast<int>(dY);

      // Make sure the interpolated point is within the raster grid (can get this kind of problem due to rounding)
      if (!bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(nRound(dXStart), nRound(dYStart), nX, nY);

      // This point is fine, so append it to the output vector
      pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // In raster-grid co-ordinates

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }

   return RTN_OK;
}
