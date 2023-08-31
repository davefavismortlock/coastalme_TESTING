/*!
 *
 * \file do_cliff_collapse.cpp
 * \brief Collapses cliffs if a critical notch depth is exceeded. Then distributes both consolidated and unconsolidated sediment from the collapse onto the shore polygons as unconsolidated talus
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
#include <assert.h>
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
   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): Calculating cliff collapse" << endl;
   
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
             dDeepWaterWavePeriod = m_VCoast[i].dGetCoastDeepWaterWavePeriod(j);
             
         // If the waves at this point are off shore, then do nothing, just move to next coast point
         if (dWaveHeightAtCoast == DBL_NODATA)
            continue;
         
         // OK we have on shore waves so calculate erosive force
         double    
             dWaveErosiveForce = pow(dWaveHeightAtCoast, WALKDEN_HALL_PARAM_1) * pow(dDeepWaterWavePeriod, WALKDEN_HALL_PARAM_2),
             dWaveEnergy = dWaveErosiveForce * m_dTimeStep * 3600;
             
         // TODO DISCUSS WITH ANDRES
         double dArbitraryConstant = 1e-5;
         dWaveEnergy *= dArbitraryConstant;

         //          assert(isfinite(dWaveEnergy));
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         // Now simulate how the coastal landform responds to this wave energy
         int nCategory = pCoastLandform->nGetLandFormCategory();
         if (nCategory == LF_CAT_CLIFF)
         {
            // This is a cliff
            CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

            // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Only consolidated sediment can have a cliff notch
            double dNotchExtension = dWaveEnergy / m_dCliffErosionResistance;

            // Extend the cliff object's erosional notch as a result of wave energy during this timestep. Note that extension may be constrained, since this-timestep extension cannot exceed the depth of sediment remaining on the cell
            pCliff->ErodeNotch(dNotchExtension);
            
            // OK, is the notch now extended enough to cause collapse (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchDepthAtCollapse))
            {
               // It is ready to collapse
               double
                   dFineCollapse = 0,
                   dSandCollapse = 0,
                   dCoarseCollapse = 0;

               // So do the cliff collapse
               nRet = nDoCliffCollapse(i, pCliff, dNotchExtension, dFineCollapse, dSandCollapse, dCoarseCollapse);
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     LogStream << m_ulIter << WARN << " problem with cliff collapse, continuing however" << endl;
               }

               // Deposit all sand and/or coarse sediment derived from this cliff collapse as unconsolidated sediment (talus)
               nRet = nDoCliffCollapseDeposition(i, pCliff,dSandCollapse, dCoarseCollapse);
               if (nRet != RTN_OK)
                  return nRet;
            }
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): total cliff collapse (m^3) = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << " (fine = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) * m_dCellArea << ", sand = " << (m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) * m_dCellArea << ", coarse = " << (m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << "), talus deposition (m^3) = " << (m_dThisIterUnconsSandCliffDeposition + m_dThisIterUnconsCoarseCliffDeposition) * m_dCellArea << " (sand = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ", coarse = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ")" << endl << endl;

   return RTN_OK;
}

/*===============================================================================================================================

Simulates cliff collapse on a single cliff object, which happens when when a notch (incised into a condsolidated sediment layer) exceeds a critical depth. This updates the cliff object, the cell 'under' the cliff object, and the polygon which contains the cliff object

===============================================================================================================================*/
int CSimulation::nDoCliffCollapse(int const nCoast, CRWCliff*pCliff, double const dNotchDeepen, double&dFineCollapse, double&dSandCollapse, double&dCoarseCollapse)
{
   // Get the cliff cell's grid coords
   int
       nX = pCliff->pPtiGetCellMarkedAsLF()->nGetX(),
       nY = pCliff->pPtiGetCellMarkedAsLF()->nGetY();
       
   // Get this cell's polygon
   int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
   if (nPoly == INT_NODATA)
   {
      // This cell isn't in a polygon
      LogStream << m_ulIter << " : in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} is not in a polygon" << endl;
      
      return RTN_ERR_CLIFF_NOT_IN_POLYGON;
   }
   
   CGeomCoastPolygon *pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
       
   // Get the elevation of the base of the notch from the cliff object
   double dNotchElev = pCliff->dGetNotchBaseElev() - m_dNotchBaseBelowSWL;

   // Get the index of the layer containing the notch (layer 0 being just above basement)
   int nNotchLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dNotchElev);
   if (nNotchLayer == ELEV_ABOVE_SEDIMENT_TOP)
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << endl;
         LogStream << m_ulIter << ": " << ERR << " cell [" << nX << "][" << nY << "] has dNotchElev (" << dNotchElev << ") above sediment top elevation (" << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << ")" << endl;
      }

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
      //      LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): all sediment removed from cliff object after cliff collapse on [" << nX << "][" << nY << "], dNotchElev = " << dNotchElev << endl;
   }

   // Now calculate the vertical depth of sediment lost in this cliff collapse. In CoastalME, all depth equivalents are assumed to be a depth upon the whole of a cell i.e. upon the area of a whole cell. The vertical depth of sediment lost in each cliff collapse is a depth upon only part of a cell, i.e. upon a fraction of a cell's area. To keep the depth of cliff collapse consistent with all other depth equivalents, weight it by the fraction of the cell's area which is being removed
   double 
      dNotchAreaFrac = dNotchDeepen / m_dCellSide,
      dAvailable = 0,
      dLost = 0,
      dFineConsLost = 0,
      dFineUnconsLost = 0,
      dSandConsLost = 0,
      dSandUnconsLost = 0,
      dCoarseConsLost = 0,
      dCoarseUnconsLost = 0;

   // LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): cell [" << nX << "][" << nY << "] before removing sediment, dGetVolEquivSedTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() << ", dGetSedimentTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << endl;

   // Now update the cell's sediment. If there are sediment layers above the notched layer, we must remove sediment from the whole depth of each layer. Again, weight the depth lost by the fraction of the cell's area which is being removed
   for (int n = nTopLayer; n > nNotchLayer; n--)
   {
      // Add in this layer's sediment, both consolidated and unconsolidated, and adjust what is left. Start with the unconsolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dFineCollapse += dLost;
         dFineUnconsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchFineLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dSandCollapse += dLost;
         dSandUnconsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchSandLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dCoarseCollapse += dLost;
         dCoarseUnconsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->IncrNotchCoarseLost(dLost);
      }

      // Now get the consolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dFineCollapse += dLost;
         dFineConsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchFineLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dSandCollapse += dLost;
         dSandConsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchSandLost(dLost);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dLost = dAvailable * dNotchAreaFrac;
         dCoarseCollapse += dLost;
         dCoarseConsLost += dLost;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->IncrNotchCoarseLost(dLost);
      }
   }

   // Now do the the consolidated layer which contains the notch. Here, we remove only part of the sediment depth
   double
       dNotchLayerTop = m_pRasterGrid->m_Cell[nX][nY].dCalcLayerElev(nNotchLayer),
       dNotchLayerThickness = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->dGetTotalThickness(),
       dNotchLayerVertFracRemoved = (dNotchLayerTop - dNotchElev) / dNotchLayerThickness;

   // Now calculate the fraction of the volume which is removed
   double dNotchLayerFracRemoved = dNotchLayerVertFracRemoved * dNotchAreaFrac;

   // Sort out the notched layer's sediment, both consolidated and unconsolidated, for this cell. First the unconsolidated sediment
   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated fine sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchFineLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated sand sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchSandLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated coarse sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->IncrNotchCoarseLost(dLost);
   }

   // Now do the consolidated sediment
   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetFine() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      // Some consolidated fine sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchFineLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetSand() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      // Some consolidated sand sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchSandLost(dLost);
   }

   dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetCoarse() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      // Some consolidated coarse sediment is available for collapse
      dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->IncrNotchCoarseLost(dLost);
   }
   
   // Update the cell's totals for cliff collapse erosion
   m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapseErosion(dFineCollapse, dSandCollapse, dCoarseCollapse);

   // Update the cell's layer elevations and d50
   m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

   // And update the cell's sea depth
   m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

   // The notch has gone
   pCliff->SetNotchDepth(0);

   // LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): cell [" << nX << "][" << nY << "] after removing sediment, dGetVolEquivSedTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() << ", dGetSedimentTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << endl << endl;
   
   // Update this-polygon totals: add to the depths of cliff collapse erosion for this polygon
   pPolygon->AddCliffCollapseErosionFine(dFineCollapse);
   pPolygon->AddCliffCollapseErosionSand(dSandCollapse);
   pPolygon->AddCliffCollapseErosionCoarse(dCoarseCollapse);
         
   // And update the this-timestep totals and the grand totals for collapse
   m_nNThisIterCliffCollapse++;
   m_nNTotCliffCollapse++;

   // Note that this gets added in to the suspended load elsewhere, so no need to do it here
   m_dThisIterCliffCollapseErosionFineUncons += dFineUnconsLost;
   m_dThisIterCliffCollapseErosionFineCons += dFineConsLost;
   
   m_dThisIterCliffCollapseErosionSandUncons += dSandUnconsLost;
   m_dThisIterCliffCollapseErosionSandCons += dSandConsLost;
   m_dThisIterCliffCollapseErosionCoarseUncons += dCoarseUnconsLost;
   m_dThisIterCliffCollapseErosionCoarseCons += dCoarseConsLost;

   return RTN_OK;
}

/*===============================================================================================================================

Redistributes the sand-sized and coarse-sized sediment from a cliff collapse onto the foreshore, as unconsolidated talus

The talus is added to the existing beach volume (i.e. to the unconsolidated sediment). The shoreline is iteratively advanced seaward until all this volume is accommodated under a Dean equilibrium profile. This equilibrium beach profile is h(y) = A * y^(2/3) where h(y) is the water depth at a distance y from the shoreline and A is a sediment-dependent scale parameter

===============================================================================================================================*/
int CSimulation::nDoCliffCollapseDeposition(int const nCoast, CRWCliff*pCliff, double const dSandCollapse, double const dCoarseCollapse)
{
   // Do we have any sand- or coarse-sized sediment to deposit?
   if ((dSandCollapse + dCoarseCollapse) < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   // OK, we have some sand- and/or coarse-sized sediment to deposit
   int
       nStartPoint = pCliff->nGetPointOnCoast(),
       nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
       
   // Get the cliff cell's grid coords
   int
       nXCliff = pCliff->pPtiGetCellMarkedAsLF()->nGetX(),
       nYCliff = pCliff->pPtiGetCellMarkedAsLF()->nGetY();    

   // Get this cell's polygon
   int nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
   if (nPoly == INT_NODATA)
   {
      // This cell isn't in a polygon
      LogStream << m_ulIter << " : in nDoCliffCollapse(), [" << nXCliff << "][" << nYCliff << "] = {" << dGridCentroidXToExtCRSX(nXCliff) << ", " << dGridCentroidYToExtCRSY(nYCliff) << "} is not in a polygon" << endl;
      
      return RTN_ERR_CLIFF_NOT_IN_POLYGON;
   }
   
   CGeomCoastPolygon *pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
       
   double
       dTotSandToDeposit = dSandCollapse,
       dTotCoarseToDeposit = dCoarseCollapse,
       dSandProp = dSandCollapse / (dSandCollapse + dCoarseCollapse),
       dCoarseProp = 1 - dSandProp;

   //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): coast = " << nCoast << ", point = " << nStartPoint << endl;

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
      nVProfileLength[n] = nRound(m_dCliffTalusMinDepositionLength);
//      dVToDepositPerProfile[n] = (dTotSandToDeposit + dTotCoarseToDeposit) / m_nCliffCollapseTalusPlanviewWidth;
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

   // Process each deposition profile
   for (int nAcross = 0; nAcross < m_nCliffCollapseTalusPlanviewWidth; nAcross++)
   {
      // Re-calculate the amount to be deposited per talus profile (this needs to change if there has been erosion of cells adding to the total to be deposited)
      for (int n = nAcross; n < m_nCliffCollapseTalusPlanviewWidth; n++)
         dVToDepositPerProfile[n] = (dTotSandToDeposit + dTotCoarseToDeposit) / (m_nCliffCollapseTalusPlanviewWidth - n);
      
      int nWidthDistSigned = nVWidthDistSigned[nAcross];

      // Get the start point of the cliff collapse deposition profile
      int nThisPoint = nStartPoint + nWidthDistSigned;

      // Is this start point valid?
      if ((nThisPoint < 0) || (nThisPoint > (nCoastSize - 1)))
      {
         // No, it is outside the grid
         if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         {                  
            LogStream << m_ulIter << ": unable to deposit sufficient unconsolidated talus from cliff collapse, since start point of Dean profile is outside the grid, nThisPoint = " << nThisPoint << endl;
            
            return RTN_ERR_CLIFFDEPOSIT;
         }
      }

      // OK, the start point of the deposition profile is inside the grid
      CGeom2DPoint
          PtStart,
          PtEnd;

      // Make the start of the deposition profile the cliff cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
      PtStart.SetX(dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetX()));
      PtStart.SetY(dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetY()));
      
      // Set the initial fraction of cliff height
      double dCliffHeightFrac = m_dMinCliffTalusHeightFrac;

      // The initial seaward offset, in cells
      int nSeawardOffset = -1;
      do
      {
         // Increase the seaward offset each time round the loop
         nSeawardOffset++;
         
         // Has the seaward offset reached the limit?
         if (nSeawardOffset > tMin(m_nXGridMax, m_nYGridMax)) 
         {
            // It has, so try again with a larger fraction of cliff height
            nSeawardOffset = 0;
            dCliffHeightFrac += 0.1;
            
            // Has the cliff height reached the limit?
            if (dCliffHeightFrac > 1)
            {
               // It has, so try again with a longer planview deposition length
               nVProfileLength[nAcross] += 10;     // TODO make this a named constant
               
               nSeawardOffset = 0;
               dCliffHeightFrac = m_dMinCliffTalusHeightFrac;
            }
               
            // Safety check: has the planview deposition length reached the limit?
            if (nVProfileLength[nAcross] > tMin(m_nXGridMax, m_nYGridMax))               
            {
               LogStream << m_ulIter << ": Unable to deposit sufficient unconsolidated talus from cliff collapse, nSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nVProfileLength[nAcross] = " << nVProfileLength[nAcross] << endl;
               
               return RTN_ERR_CLIFFDEPOSIT;
            }
         }            

         // Now construct a deposition collapse profile from the start point, it is one cell longer than the specified length because it includes the cliff point in the profile. Calculate its length in external CRS units, the way it is done here is approximate but probably OK
         double dThisProfileLength = (nVProfileLength[nAcross] + nSeawardOffset + 1) * m_dCellSide;

         // Get the end point of this coastline-normal line
         CGeom2DIPoint PtiEnd; // In grid CRS
         int nRtn = nGetCoastNormalEndPoint(nCoast, nThisPoint, nCoastSize, &PtStart, dThisProfileLength, &PtEnd, &PtiEnd);
         if (nRtn == RTN_ERR_PROFILE_ENDPOINT_IS_OFFGRID)
         {
            LogStream << m_ulIter << ": Unable to deposit sufficient unconsolidated talus from cliff collapse, since end point of Dean profile has reached the edge of the grid with nSeawardOffset = " << nSeawardOffset << endl;
            
            return nRtn;
         }

         if (nRtn == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
         {
            LogStream << m_ulIter << ": Unable to deposit sufficient unconsolidated talus from cliff collapse, could not find a solution for the end point of the Dean profile" << endl;
            
            return nRtn;
         }

         // OK, both the start and end points of this depoition profile are within the grid
         //         LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): nWidthDistSigned = " << nWidthDistSigned << " cliff collapse profile from " << PtStart.dGetX() << ", " << PtStart.dGetY() << " to " << PtEnd.dGetX() << ", " << PtEnd.dGetY() << " with length (inc. cliff point) = " << dThisProfileLength << endl;

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
         
         // Calculate the existing elevation for all points along the deposition profile
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

         // Now calculate the elevation of the talus top at the shoreline
         double
             dCliffTopElev = dVProfileNow[0],
             dCliffBaseElev = dVProfileNow[1],
             dCliffHeight = dCliffTopElev - dCliffBaseElev,
             dTalusTopElev = dCliffBaseElev + (dCliffHeight * dCliffHeightFrac);

         //         LogStream << "Elevations: cliff top = " << dCliffTopElev << " cliff base = " << dCliffBaseElev << " talus top = " << dTalusTopElev << endl;

         //          if (dCliffTopElev < dCliffBaseElev)
         //             LogStream << "*** ERROR, cliff top is lower than cliff base" << endl;

         // Next calculate the talus slope length in external CRS units, this is approximate but probably OK
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
         if (dTotElevDiff >= dVToDepositPerProfile[nAcross])
         {
            // No it doesn't, so try again with a larger seaward offset
            // LogStream << m_ulIter << ": nSeawardOffset = " << nSeawardOffset << " dTotElevDiff = " << dTotElevDiff << " dVToDepositPerProfile[nAcross] = " << dVToDepositPerProfile[nAcross] << endl;
            
            continue;
         }

         // Yes it does
         //          LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): cliff collapse at [" << VCellsUnderProfile[0].nGetX() << "][" << VCellsUnderProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(VCellsUnderProfile[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(VCellsUnderProfile[0].nGetY()) << "} offset SUFFICIENT with nSeawardOffset = " << nSeawardOffset << endl;
         //          LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): dTotElevDiff = " << dTotElevDiff << " dVToDepositPerProfile[nAcross] = " << dVToDepositPerProfile[nAcross] << endl;

//         double dPropToDeposit = dVToDepositPerProfile[nAcross] / dTotElevDiff;
         //         LogStream << "dPropToDeposit = " << dPropToDeposit << endl;

         //          double
         //             dDepositedCheck = 0,
         //             dRemovedCheck = 0;

         // Process all cells in this profile
         for (int n = 0; n < nRasterProfileLength; n++)
         {
            // Have we deposited all that we need to?
            if ((dTotSandToDeposit + dTotCoarseToDeposit) <= 0)
            {
               LogStream << m_ulIter << ": (dTotSandToDeposit + dTotCoarseToDeposit) = " << (dTotSandToDeposit + dTotCoarseToDeposit)<< endl;
               break;
            }
            
            // OK, we still have some talus left to deposit
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
                  dPotentialSandToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dSandProp;         //  * dPropToDeposit;
                  dPotentialSandToDeposit = tMin(dPotentialSandToDeposit, dTotSandToDeposit);

                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSand(dPotentialSandToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotSandToDeposit -= dPotentialSandToDeposit;
                  //                  dDepositedCheck += dPotentialSandToDeposit;
               }

               double dPotentialCoarseToDeposit = 0;
               if (dTotCoarseToDeposit > 0)
               {
                  dPotentialCoarseToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dCoarseProp;        //  * dPropToDeposit;
                  dPotentialCoarseToDeposit = tMin(dPotentialCoarseToDeposit, dTotCoarseToDeposit);

                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarse(dPotentialCoarseToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dTotCoarseToDeposit -= dPotentialCoarseToDeposit;
                  //                  dDepositedCheck += dPotentialCoarseToDeposit;
               }

               // Now update the cell's layer elevations
               m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

               // Update the cell's sea depth
               m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

               // Update the cell's talus deposition, and total talus deposition, values
               m_pRasterGrid->m_Cell[nX][nY].AddSandTalusDeposition(dPotentialSandToDeposit);
               m_pRasterGrid->m_Cell[nX][nY].AddCoarseTalusDeposition(dPotentialCoarseToDeposit);

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

                  // And increment the per-timestep total for fine sediment eroded during cliff collapse deposition (note that this gets added in to the suspended load elsewhere, so no need to do it here)
                  m_dThisIterCliffCollapseFineErodedDuringDeposition += dFine;
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

                  // And increment the per-timestep total for sand sediment eroded during cliff collapse deposition
                  m_dThisIterCliffCollapseSandErodedDuringDeposition += dSand;
                  
                  // Add this to the total of sand-sized talus to be deposited
                  dTotSandToDeposit += dSand;
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

                  // And increment the per-timestep total for coarse sediment eroded during cliff collapse deposition
                  m_dThisIterCliffCollapseCoarseErodedDuringDeposition += dCoarse;
                  
                  // Add this to the total of coarse talus to be deposited
                  dTotCoarseToDeposit += dCoarse;                  
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
         }     // All cells in this profile
         //        LogStream << endl;
         //        LogStream << "Profile done, dDepositedCheck = " << dDepositedCheck << " dRemovedCheck = " << dRemovedCheck << endl;
         
         LogStream << m_ulIter << ": XXXX" << endl;
         break;
         
      } while (true);      // The seaward offset loop
   }     // Process each deposition profile

   // Safety check for sand sediment
   if (! bFPIsEqual(dTotSandToDeposit, 0.0, TOLERANCE))
   {
      LogStream << ERR << m_ulIter << ": dTotSandToDeposit = " << dTotSandToDeposit << " SET TO ZERO" << endl;
      dTotSandToDeposit = 0;
   }

   // Ditto for coarse sediment
   if (! bFPIsEqual(dTotCoarseToDeposit, 0.0, TOLERANCE))
   {
      LogStream << ERR << m_ulIter << ": dTotCoarseToDeposit = " << dTotCoarseToDeposit << " SET TO ZERO" << endl;
      dTotCoarseToDeposit = 0;
   }

   // Store the total depths of cliff collapse deposition for this polygon
   pPolygon->AddCliffCollapseUnconsSandDeposition(dSandCollapse);
   pPolygon->AddCliffCollapseUnconsCoarseDeposition(dCoarseCollapse);
   
   // Increment this-timestep totals for cliff collapse deposition
   m_dThisIterUnconsSandCliffDeposition += dSandCollapse;
   m_dThisIterUnconsCoarseCliffDeposition += dCoarseCollapse;

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
      if (! bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(nRound(dXStart), nRound(dYStart), nX, nY);

      // This point is fine, so append it to the output vector
      pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // In raster-grid co-ordinates

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }

   return RTN_OK;
}
