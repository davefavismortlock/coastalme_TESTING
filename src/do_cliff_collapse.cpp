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

//===============================================================================================================================
//! Update accumulated wave energy in coastal landform objects
//===============================================================================================================================
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

         // First get wave energy for the coastal landform object
         double dWaveHeightAtCoast = m_VCoast[i].dGetCoastWaveHeight(j);
             
         // If the waves at this point are off-shore, then do nothing, just move to next coast point
         if (bFPIsEqual(dWaveHeightAtCoast, DBL_NODATA, TOLERANCE))
            continue;
         
         // OK we have on-shore waves so get the previously-calculated wave energy
         double dWaveEnergy = m_VCoast[i].dGetWaveEnergyAtBreaking(j);    
//          assert(isfinite(dWaveEnergy));

         // And save the accumulated value
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         // Now simulate how the coastal landform responds to this wave energy
         int nCategory = pCoastLandform->nGetLandFormCategory();
         if (nCategory == LF_CAT_CLIFF)
         {
            // This is a cliff
            CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

            // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Only consolidated sediment can have a cliff notch
            double dNotchExtension = dWaveEnergy / m_dCliffErosionResistance;

            // Deepen the cliff object's erosional notch as a result of wave energy during this timestep. Note that notch deepening may be constrained, since this-timestep notch extension cannot exceed the length (i.e. cellside minus notch depth) of sediment remaining on the cell
            pCliff->DeepenErosionalNotch(dNotchExtension);
            
            // OK, is the notch now extended enough to cause collapse (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchDepthAtCollapse))
            {
               // It is ready to collapse
               double
                  dCliffElevPreCollapse = 0,
                  dCliffElevPostCollapse = 0,
                  dFineCollapse = 0,
                  dSandCollapse = 0,
                  dCoarseCollapse = 0;
                  
               // So do the cliff collapse
               nRet = nDoCliffCollapse(i, pCliff, dFineCollapse, dSandCollapse, dCoarseCollapse, dCliffElevPreCollapse, dCliffElevPostCollapse);
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     LogStream << m_ulIter << WARN << " problem with cliff collapse, continuing however" << endl;
               }

               // Deposit all sand and/or coarse sediment derived from this cliff collapse as unconsolidated sediment (talus)
               nRet = nDoCliffCollapseDeposition(i, pCliff,dSandCollapse, dCoarseCollapse, dCliffElevPreCollapse, dCliffElevPostCollapse);
               if (nRet != RTN_OK)
                  return nRet;
            }
         }
      }
   }

   // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
   //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): total cliff collapse (m^3) = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << " (fine = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) * m_dCellArea << ", sand = " << (m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) * m_dCellArea << ", coarse = " << (m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << "), talus deposition (m^3) = " << (m_dThisIterUnconsSandCliffDeposition + m_dThisIterUnconsCoarseCliffDeposition) * m_dCellArea << " (sand = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ", coarse = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ")" << endl << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Simulates cliff collapse on a single cliff object, which happens when when a notch (incised into a condsolidated sediment layer) exceeds a critical depth. This updates the cliff object, the cell 'under' the cliff object, and the polygon which contains the cliff object
//===============================================================================================================================
int CSimulation::nDoCliffCollapse(int const nCoast, CRWCliff*pCliff, double& dFineCollapse, double& dSandCollapse, double& dCoarseCollapse, double& dPreCollapseCliffElev, double& dPostCollapseCliffElev)
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
      // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      //    LogStream << m_ulIter << ": " << ERR << " cell [" << nX << "][" << nY << "] has dNotchElev (" << dNotchElev << ") above sediment top elevation (" << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << ")" << endl;

      return RTN_ERR_CLIFFNOTCH;
   }

   // Flag the coastline cliff object as having collapsed
   pCliff->SetCliffCollapsed();

   int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

   // Safety check
   if (nTopLayer == INT_NODATA)
      return RTN_ERR_NO_TOP_LAYER;

   // Set the base of the collapse (see above)
   pCliff->SetNotchBaseElev(dNotchElev);
   
   // Set flags to say that the top layer has changed
   m_bConsChangedThisIter[nTopLayer] = true;
   m_bUnconsChangedThisIter[nTopLayer] = true;
   
   // Get the pre-collapse cliff elevation
   dPreCollapseCliffElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

   // Now calculate the vertical depth of sediment lost in this cliff collapse. In CoastalME, all depth equivalents are assumed to be a depth upon the whole of a cell i.e. upon the area of a whole cell. So to keep the depth of cliff collapse consistent with all other depth equivalents, weight it by the fraction of the cell's area which is being removed
   double 
      dAvailable = 0,
      dFineConsLost = 0,
      dFineUnconsLost = 0,
      dSandConsLost = 0,
      dSandUnconsLost = 0,
      dCoarseConsLost = 0,
      dCoarseUnconsLost = 0;

   // Now update the cell's sediment. If there are sediment layers above the notched layer, we must remove sediment from the whole depth of each layer
   for (int n = nTopLayer; n > nNotchLayer; n--)
   {
      // Start with the unconsolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);         
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);         
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);         
      }

      // Now get the consolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchFineLost();
      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchFineLost(0);         
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchSandLost();
      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchSandLost(0);         
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);         
      }
   }

   // Now sort out the sediment lost from the consolidated layer into which the erosional notch was incised
   double
       dNotchLayerTop = m_pRasterGrid->m_Cell[nX][nY].dCalcLayerElev(nNotchLayer),
       dNotchLayerThickness = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->dGetTotalThickness(),
       dNotchLayerFracRemoved = (dNotchLayerTop - dNotchElev) / dNotchLayerThickness;

   // Sort out the notched layer's sediment, both consolidated and unconsolidated, for this cell. First the unconsolidated sediment
   double dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated fine sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);
   }

   double dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated sand sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);
   }

   double dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      // Some unconsolidated coarse sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);
   }

   // Now do the consolidated sediment
   dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchFineLost();
   if (dAvailable > 0)
   {
      // Some consolidated fine sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchFineLost(0);
   }

   dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchSandLost();
   if (dAvailable > 0)
   {
      // Some consolidated sand sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchSandLost(0);
   }

   dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchCoarseLost();
   if (dAvailable > 0)
   {
      // Some consolidated coarse sediment is available for collapse
      double dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);
   }
   
   // Update the cell's totals for cliff collapse erosion
   m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapseErosion(dFineCollapse, dSandCollapse, dCoarseCollapse);

   // Update the cell's layer elevations and d50
   m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
   
   // Get the post-collapse cliff elevation
   dPostCollapseCliffElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

   // And update the cell's sea depth
   m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

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

//===============================================================================================================================
//! Redistributes the sand-sized and coarse-sized sediment from a cliff collapse onto the foreshore, as unconsolidated talus. The talus is added to the existing beach volume (i.e. to the unconsolidated sediment). The shoreline is iteratively advanced seaward until all this volume is accommodated under a Dean equilibrium profile. This equilibrium beach profile is h(y) = A * y^(2/3) where h(y) is the water depth at a distance y from the shoreline and A is a sediment-dependent scale parameter
//===============================================================================================================================
int CSimulation::nDoCliffCollapseDeposition(int const nCoast, CRWCliff*pCliff, double const dSandFromCollapse, double const dCoarseFromCollapse, double const dPreCollapseCliffElev, double const dPostCollapseCliffElev)
{
   // Check: is there some sand- or coarse-sized sediment to deposit?
   if ((dSandFromCollapse + dCoarseFromCollapse) < SEDIMENT_ELEV_TOLERANCE)
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
       
   // Use this to make sure that, if we have both sand and coarse to deposit, we don't drop all the sand on a cell and then be unable to deposit any coarse
   double dSandProp = dSandFromCollapse / (dSandFromCollapse + dCoarseFromCollapse);

   // OK, now set up the planview sequence talus deposition. First we deposit to create a Dean profile starting from the cliff collapse cell. Then we deposit along two profiles which start from the coast cells on either side of the cliff collapse cell, and then on two profiles starting on the next two "outside" coast cells, etc. However, we do have a "preferred" talus width
   int nTalusWidth = nConvertMetresToNumCells(m_dCliffDepositionPlanviewWidth);
   
   // The default talus collapse must be an odd number of cells in width i.e. centred on the cliff collpase cell (but only if we are not at the end of the coast)
   if ((nTalusWidth % 2) == 0)
      nTalusWidth++;
   
   // Calculate the amount to be deposited on each talus profile, assuming the "preferred" talus width
   double dProfileSandToDeposit = dSandFromCollapse / nTalusWidth;
   double dProfileCoarseToDeposit = dCoarseFromCollapse / nTalusWidth;
   
   // This holds the valid coast start points for each Dean profile
   vector<int> VnTalusProfileCoastStartPoint(nCoastSize);
   
   // This is the coast start point for the first Dean profile
   VnTalusProfileCoastStartPoint[0] = nStartPoint;
   int 
      nn = 1,
      nSigned = 1,
      nCount = 1;
   do
   {
      int nTmpPoint;
      
      if ((nCount % 2) != 0)
         nTmpPoint = nStartPoint + nSigned;
      else
      {
         nTmpPoint = nStartPoint - nSigned;
         nSigned++;
      }
      
      // Is this start point valid?
      if ((nTmpPoint < 0) || (nTmpPoint > (nCoastSize - 1)))
      {
         // No, it is outside the grid, so find another 
         nCount++;         
         continue;
      }
      else
      {
         // It is valid
         VnTalusProfileCoastStartPoint[nn] = nTmpPoint;
         nn++;
         nCount++;
      }      
   }
   while (nn < static_cast<int>(VnTalusProfileCoastStartPoint.size()));

   bool 
      bHitFirstCoastPoint = false,
      bHitLastCoastPoint = false;
   double
      dTotSandToDeposit = dSandFromCollapse,
      dTotCoarseToDeposit = dCoarseFromCollapse,
      dTotSandDeposited = 0,
      dTotCoarseDeposited = 0;
      
   // Process each deposition profile
   for (int nAcross = 0; nAcross < static_cast<int>(VnTalusProfileCoastStartPoint.size()); nAcross++)
   {
      // This holds the minimum planview length of the Dean profile (the initial length will be increased later if we can't deposit sufficient talus)
      int nTalusProfileLen = nConvertMetresToNumCells(m_dCliffTalusMinDepositionLength);
   
      double
         dProfileSandDeposited = 0,
         dProfileCoarseDeposited = 0;
         
      bool
         bSandDeposition = false,
         bSandDone = false,
         bCoarseDeposition = false,
         bCoarseDone = false;
      
      if (dTotSandToDeposit > 0)
      {
         bSandDeposition = true;
         
         if (bFPIsEqual(dTotSandToDeposit - dTotSandDeposited, 0.0, MASS_BALANCE_TOLERANCE))
            bSandDone = true;         
      }
      else
         bSandDone = true;
         
      if (dTotCoarseToDeposit > 0)
      {
         bCoarseDeposition = true;
         
         if (bFPIsEqual(dTotCoarseToDeposit - dTotCoarseDeposited, 0.0, MASS_BALANCE_TOLERANCE))
            bCoarseDone = true;
      }
      else
         bCoarseDone = true;

      if (bSandDone && bCoarseDone)
      {
         // LogStream << m_ulIter << ": break 2 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << endl << "\tnAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl;
         
         break;
      }
         
      // Get the start point of this cliff collapse deposition profile
      int nThisPoint = VnTalusProfileCoastStartPoint[nAcross];
      
      // Are we at the start or end of the coast?
      if (nThisPoint == 0)
         bHitFirstCoastPoint = true;
      
      if (nThisPoint == nCoastSize-1)
         bHitLastCoastPoint = true;

      CGeom2DPoint
          PtStart,
          PtEnd;

      // Make the start of the deposition profile the cliff cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
      PtStart.SetX(dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetX()));
      PtStart.SetY(dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetY()));
      
      // Set the initial fraction of cliff height
      double dCliffHeightFrac = m_dMinCliffTalusHeightFrac;
      
      bool
         bSandDoneOnProfile = false,
         bSandDepositionOnProfile = false,
         bCoarseDoneOnProfile = false,
         bCoarseDepositionOnProfile = false;

      // The initial seaward offset, in cells
      int nSeawardOffset = 0;
      do
      {
         if (bSandDoneOnProfile && bCoarseDoneOnProfile)   
         {
            // LogStream << m_ulIter << ": break 3 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl;
            
            break;
         }

         // Need to deposit more on this profile: increase the seaward offset each time round the loop
         nSeawardOffset++;
         
         // Has the seaward offset reached the (rather arbitrary) limit?
         if (nSeawardOffset == 20)           // tMin(m_nXGridMax, m_nYGridMax)) 
         {
            // It has, so try again with a larger fraction of cliff height
            nSeawardOffset = 1;
            dCliffHeightFrac += CLIFF_COLLAPSE_HEIGHT_INCREMENT;
            
            // Has the cliff height reached the limit?
            if (bFPIsEqual(dCliffHeightFrac, 1.0, TOLERANCE))
            {
               // It has, so try again with a longer planview deposition length
               nTalusProfileLen += CLIFF_COLLAPSE_LENGTH_INCREMENT;
               
               nSeawardOffset = 1;
//               dCliffHeightFrac = m_dMinCliffTalusHeightFrac;
            }
         }
         
         if (bHitFirstCoastPoint && bHitLastCoastPoint)
         {
            // Uh-oh, we've reached both ends of the coast (!) and we can't increase anything any more
            LogStream << m_ulIter << ": unable to deposit sufficient unconsolidated talus from cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl;
            
            return RTN_ERR_CLIFFDEPOSIT;
         }

         // Now construct a deposition collapse profile from the start point, it is one cell longer than the specified length because it includes the cliff point in the profile. Calculate its length in external CRS units, the way it is done here is approximate but probably OK
         double dThisProfileLength = (nTalusProfileLen + nSeawardOffset + 1) * m_dCellSide;

         // Get the end point of this coastline-normal line
         CGeom2DIPoint PtiEnd;            // In grid CRS
         int nRtn = nGetCoastNormalEndPoint(nCoast, nThisPoint, nCoastSize, &PtStart, dThisProfileLength, &PtEnd, &PtiEnd);
         if (nRtn == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
         {
            LogStream << m_ulIter << ": unable to deposit sufficient unconsolidated talus from cliff collapse, could not find a solution for the end point of the Dean profile" << endl;
            
            return nRtn;
         }
         
         // OK, both the start and end points of this deposition profile are within the grid
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
             dCliffHeight = dPreCollapseCliffElev - dPostCollapseCliffElev,
             dTalusTopElev = dPostCollapseCliffElev + (dCliffHeight * dCliffHeightFrac);

         //         LogStream << "Elevations: cliff top = " << dPreCollapseCliffElev << " cliff base = " << dPostCollapseCliffElev << " talus top = " << dTalusTopElev << endl;

         //          if (dPreCollapseCliffElev < dPostCollapseCliffElev)
         //             LogStream << "*** ERROR, cliff top is lower than cliff base" << endl;

         // Next calculate the talus slope length in external CRS units, this is approximate but probably OK
         double dTalusSlopeLength = dThisProfileLength - ((nSeawardOffset - 1) * m_dCellSide);

         // If user has not supplied a value for m_dCliffDepositionA, then solve for dA so that the elevations at end of the existing profile, and at the end of the Dean equilibrium profile, are the same
         double dA = 0;
         if (bFPIsEqual(m_dCliffDepositionA, 0.0, TOLERANCE))
            dA = m_dCliffDepositionA;
         else
            dA = (dTalusTopElev - dVProfileNow[nRasterProfileLength - 1]) / pow(dTalusSlopeLength, DEAN_POWER);

         double dInc = dTalusSlopeLength / (nRasterProfileLength - nSeawardOffset - 2);
         vector<double> dVDeanProfile(nRasterProfileLength);

         // Calculate the Dean equilibrium profile of the talus h(y) = A * y^(2/3) where h(y) is the distance below the talus-top elevation (the highest point in the Dean profile) at a distance y from the cliff (the landward start of the profile)
         CalcDeanProfile(&dVDeanProfile, dInc, dTalusTopElev, dA, true, nSeawardOffset, dTalusTopElev);

         // Get the total difference in elevation between the two profiles (Dean profile - present profile). Since we want the Dean profile to be higher than the present profile, a good result is a +ve number
         double dTotElevDiff = dSubtractProfiles(&dVDeanProfile, &dVProfileNow, &bVProfileValid);

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
         if (dTotElevDiff < (dProfileSandToDeposit + dProfileCoarseToDeposit))
         {
            // No it doesn't, so try again with a larger seaward offset and/or a longer Dean profile length
            // LogStream << m_ulIter << ": nSeawardOffset = " << nSeawardOffset << " dTotElevDiff = " << dTotElevDiff << " VdTalusToDepositPerProfile[nAcross] = " << VdTalusToDepositPerProfile[nAcross] << endl;
            
            continue;
         }
         
         // Yes it does, so set the switches for sand and coarse deposition
         if (! bFPIsEqual(dProfileSandToDeposit, 0.0, MASS_BALANCE_TOLERANCE))
         {
            bSandDepositionOnProfile = true;
            bSandDoneOnProfile = false;
         }
         else
         {
            bSandDepositionOnProfile = false;
            bSandDoneOnProfile = true;
         }
         
         if (! bFPIsEqual(dProfileCoarseToDeposit, 0.0, MASS_BALANCE_TOLERANCE))
         {
            bCoarseDepositionOnProfile = true;
            bCoarseDoneOnProfile = false;
         }
         else
         {
            bCoarseDepositionOnProfile = false;
            bCoarseDoneOnProfile = true;
         }
         
         // OK, now process all cells in this profile, including the first one (which is where the cliff collapse occurred)
         for (int n = 0; n < nRasterProfileLength; n++)
         {
            // Are we depositing sand talus sediment on this profile?
            if (bSandDepositionOnProfile)
            {
               if (bFPIsEqual(dProfileSandToDeposit - dProfileSandDeposited, 0.0, MASS_BALANCE_TOLERANCE))
                  bSandDoneOnProfile = true;
               else
                  bSandDoneOnProfile = false;
            }

            // Are we depositing coarse talus sediment on this profile?
            if (bCoarseDepositionOnProfile)
            {
               if (bFPIsEqual(dProfileCoarseToDeposit - dProfileCoarseDeposited, 0.0, MASS_BALANCE_TOLERANCE))
                  bCoarseDoneOnProfile = true;
               else
                  bCoarseDoneOnProfile = false;
            }

            // If we have deposited enough, then break out of the loop
            if (bSandDoneOnProfile && bCoarseDoneOnProfile)   
            {
               // LogStream << m_ulIter << ": break 1 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl;
               
               break;
            }

            // Nope, we still have some talus left to deposit
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

            // Only do deposition on this cell if its elevation is below the Dean elevation
            if (dVDeanProfile[n] > dVProfileNow[n])
            {
              // At this point along the profile, the Dean profile is higher than the present profile. So we can deposit some sediment on this cell
               double dSandToDeposit = 0;
               if (bSandDepositionOnProfile)
               {
                  dSandToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dSandProp;
                  dSandToDeposit = tMin(dSandToDeposit, (dProfileSandToDeposit - dProfileSandDeposited), (dTotSandToDeposit - dTotSandDeposited));

                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;
                  
                  dProfileSandDeposited += dSandToDeposit;
                  dTotSandDeposited += dSandToDeposit;
               }

               double dCoarseToDeposit = 0;
               if (bCoarseDepositionOnProfile)
               {
                  dCoarseToDeposit = (dVDeanProfile[n] - dVProfileNow[n]);
                  dCoarseToDeposit = tMin(dCoarseToDeposit, (dProfileCoarseToDeposit - dProfileCoarseDeposited), (dTotCoarseToDeposit - dTotCoarseDeposited));

                  m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseToDeposit);

                  // Set the changed-this-timestep switch
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  dProfileCoarseDeposited += dCoarseToDeposit;
                  dTotCoarseDeposited += dCoarseToDeposit;
               }

               // Now update the cell's layer elevations
               m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

               // Update the cell's sea depth
               m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

               // Update the cell's talus deposition, and total talus deposition, values
               m_pRasterGrid->m_Cell[nX][nY].AddSandTalusDeposition(dSandToDeposit);
               m_pRasterGrid->m_Cell[nX][nY].AddCoarseTalusDeposition(dCoarseToDeposit);

               // And set the landform category
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_TALUS);
            }

//             else if (dVDeanProfile[n] < dVProfileNow[n])
//             {
//                // Here, the Dean profile is lower than the present profile, so we must remove some sediment from this cell
//                double dThisLowering = dVProfileNow[n] - dVDeanProfile[n];
// 
//                // Find out how much sediment we have available on this cell
//                double dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
//                double dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
//                double dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
// 
//                // Now partition the total lowering for this cell between the three size fractions: do this by relative erodibility
//                int nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0);
//                int nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0);
//                int nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);
// 
//                double dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);
// 
//                if (nFineWeight)
//                {
//                   // Erode some fine-sized sediment
//                   double dFineLowering = (m_dFineErodibilityNormalized * dThisLowering) / dTotErodibility;
// 
//                   // Make sure we don't get -ve amounts left on the cell
//                   double dFine = tMin(dExistingAvailableFine, dFineLowering);
//                   double dRemaining = dExistingAvailableFine - dFine;
// 
//                   // Set the value for this layer
//                   m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dRemaining);
// 
//                   // And set the changed-this-timestep switch
//                   m_bUnconsChangedThisIter[nTopLayer] = true;
// 
//                   // And increment the per-timestep total for fine sediment eroded during cliff collapse deposition (note that this gets added in to the suspended load elsewhere, so no need to do it here)
//                   m_dThisIterCliffCollapseFineErodedDuringDeposition += dFine;
//                }
// 
//                if (nSandWeight)
//                {
//                   // Erode some sand-sized sediment
//                   double dSandLowering = (m_dSandErodibilityNormalized * dThisLowering) / dTotErodibility;
// 
//                   // Make sure we don't get -ve amounts left on the source cell
//                   double dSandToErode = tMin(dExistingAvailableSand, dSandLowering);
//                   double dRemaining = dExistingAvailableSand - dSandToErode;
// 
//                   // Set the value for this layer
//                   m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dRemaining);
// 
//                   // Set the changed-this-timestep switch
//                   m_bUnconsChangedThisIter[nTopLayer] = true;
// 
//                   // And increment the per-timestep total for sand sediment eroded during cliff collapse deposition
//                   m_dThisIterCliffCollapseSandErodedDuringDeposition += dSandToErode;
//                   
//                   // Increase the sand deposition target
//                   dTotSandToDeposit += dSandToErode;
//                   
//                   // And set switches
//                   bSandDeposition = true;
//                   bSandDepositionOnProfile = true;
//                   
//                }
// 
//                if (nCoarseWeight)
//                {
//                   // Erode some coarse-sized sediment
//                   double dCoarseLowering = (m_dCoarseErodibilityNormalized * dThisLowering) / dTotErodibility;
// 
//                   // Make sure we don't get -ve amounts left on the source cell
//                   double dCoarseToErode = tMin(dExistingAvailableCoarse, dCoarseLowering);
//                   double dRemaining = dExistingAvailableCoarse - dCoarseToErode;
// 
//                   // Set the value for this layer
//                   m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dRemaining);
// 
//                   // Set the changed-this-timestep switch
//                   m_bUnconsChangedThisIter[nTopLayer] = true;
// 
//                   // And increment the per-timestep total for coarse sediment eroded during cliff collapse deposition
//                   m_dThisIterCliffCollapseCoarseErodedDuringDeposition += dCoarseToErode;
//                   
//                   // Increase the coarse deposition target
//                   dTotCoarseToDeposit += dCoarseToErode;
//                   
//                   // And set switches
//                   bCoarseDeposition = true;
//                   bCoarseDepositionOnProfile = true;
//                   
//                }
// 
//                // Recalculate the elevation of every layer
//                m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
// 
//                // And update the cell's sea depth
//                m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
//             }
         }     // All cells in this profile 
         
         // OK we have either processed all cells in this profile, or we have deposited enough talus sediment on this profile
         if (bSandDeposition)
         {
            double dNotDeposited = dProfileSandToDeposit - dProfileSandDeposited;
            if (! bFPIsEqual(dNotDeposited, 0.0, MASS_BALANCE_TOLERANCE))
               dTotSandToDeposit += dNotDeposited;
         }
            
         if (bCoarseDeposition)
         {
            double dNotDeposited = dProfileCoarseToDeposit - dProfileCoarseDeposited;
            if (! bFPIsEqual(dNotDeposited, 0.0, MASS_BALANCE_TOLERANCE))
               dTotCoarseToDeposit += dNotDeposited;
         }
         
         break;         
      } 
      while (true);      // The seaward offset etc. loop
      
      // LogStream << m_ulIter << ": after while-seaward offset loop for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl;
      
      // Have we deposited enough for this cliff collapse?
      if (bSandDeposition)
      {
         if (bFPIsEqual(dTotSandDeposited - dTotSandToDeposit, 0.0, MASS_BALANCE_TOLERANCE))
            bSandDone = true;
      }
      
      if (bCoarseDeposition)
      {
         if (bFPIsEqual(dTotCoarseDeposited - dTotCoarseToDeposit, 0.0, MASS_BALANCE_TOLERANCE))
            bCoarseDone = true;
      }
      
      if (bSandDone && bCoarseDone)
      {
         // LogStream << m_ulIter << ": break after all cells in profile loop for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFrac = " << dCliffHeightFrac << " nAcross = " << nAcross << endl << "\tdProfileSandToDeposit = " << dProfileSandToDeposit << " dProfileSandDeposited = " << dProfileSandDeposited << " dTotSandToDeposit = " << dTotSandToDeposit << " dTotSandDeposited = " << dTotSandDeposited << endl << "\tdProfileCoarseToDeposit = " << dProfileCoarseToDeposit << " dProfileCoarseDeposited = " << dProfileCoarseDeposited << " dTotCoarseToDeposit = " <<  dTotCoarseToDeposit << " dTotCoarseDeposited = " << dTotCoarseDeposited << endl << endl;
      
         break;
      }      
   }     // Process each deposition profile

   // Safety check for sand sediment
   if (! bFPIsEqual((dTotSandToDeposit - dTotSandDeposited), 0.0, MASS_BALANCE_TOLERANCE))
      LogStream << ERR << m_ulIter << ": non-zero dTotSandToDeposit = " << dTotSandToDeposit << endl;

   // Ditto for coarse sediment
   if (! bFPIsEqual((dTotCoarseToDeposit - dTotCoarseDeposited), 0.0, MASS_BALANCE_TOLERANCE))
      LogStream << ERR << m_ulIter << ": non-zero dTotCoarseToDeposit = " << dTotCoarseToDeposit << endl;

   // Store the total depths of cliff collapse deposition for this polygon
   pPolygon->AddCliffCollapseUnconsSandDeposition(dTotSandDeposited);
   pPolygon->AddCliffCollapseUnconsCoarseDeposition(dTotCoarseDeposited);
   
   // Increment this-timestep totals for cliff collapse deposition
   m_dThisIterUnconsSandCliffDeposition += dTotSandDeposited;
   m_dThisIterUnconsCoarseCliffDeposition += dTotCoarseDeposited;

   return RTN_OK;
}

//===============================================================================================================================
//! Given the start and end points of a cliff-collapse normal profile, returns an output vector of cells which are 'under' the vector line
//===============================================================================================================================
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
