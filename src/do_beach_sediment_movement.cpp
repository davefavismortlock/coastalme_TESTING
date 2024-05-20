/*!
 *
 * \file do_beach_sediment_movement.cpp
 * \brief Does between-polygon actual (supply-limited) redistribution of transported beach sediment
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
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
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::to_string;

#include <algorithm>
using std::stable_sort;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

//===============================================================================================================================

//! Function used to sort polygons before doing the polygon-to-polygon source-target pattern. For both LH and RH arguments, the first value is the polygon coast ID, the second value is the down- or up-coast direction, and subsequent numbers are adjacent polygon coastIDs in that direction. If the first argument must be ordered before the second, return true

//===============================================================================================================================
bool bPolygonAndAdjCompare(const vector<int>& nVLeft, const vector<int>& nVRight)
{
   // For safety, check that the LHS polygon has at least one adjacent polygon (it should have, apart from the bad situation where just one big polygon is created)
   if ((nVLeft.size() >= 3) && (nVRight.size() >= 3))
   {
      // Polygons at the grid edge are processed last, so put LHS grid-edge polygons on the RHS
      if (nVLeft[2] == INT_NODATA)
         return false;

      // Polygons at the grid edge are processed last, so keep RHS grid-edge polygons where they are
      if (nVRight[2] == INT_NODATA)
         return true;

      // Now sort out polygon-to-polygon dependencies. We need to put 'target' polygons after 'source' polygons, so that the source is processed before the target. So does the LHS polygon have the RHS polygon as one of its adjacent polygons?
      for (unsigned int n = 2; n < nVLeft.size(); n++)
      {
         if (nVRight[0] == nVLeft[n])
            // It does, so keep the existing sequence
            return true;
      }

      // Does the RHS polygon have the LHS polygon as one of its adjacent polygons?
      for (unsigned int n = 2; n < nVRight.size(); n++)
      {
         if (nVLeft[0] == nVRight[n])
            // It does, so swap them
            return false;
      }
   }

   bool bDownCoast = nVLeft[1];
   if (bDownCoast)
      // Sediment going down-coast
      return nVLeft < nVRight;
   else
      // Sediment going up-coast
      return nVLeft > nVRight;

   // Default return value, should never get here
   return true;
}

//===============================================================================================================================
//! Does between-polygon and within-polygon actual (supply-limited) redistribution of transported beach sediment
//===============================================================================================================================
int CSimulation::nDoAllActualBeachErosionAndDeposition(void)
{
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nRet;

      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): Calculating unconsolidated sediment transport" << endl;
      
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         WritePolygonShareTable(nCoast);         
         WritePolygonPreExistingSediment(nCoast);         
         WritePolygonShorePlatformErosion(nCoast);         
         WritePolygonCliffCollapseErosion(nCoast);
      }      

      // Update the values of pre-existing unconsolidated sediment, for all three size classes, to include unconsolidated sediment derived from platform erosion and/or cliff collapse
      AllPolygonsUpdateStoredUncons(nCoast);

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         WritePolygonSedimentBeforeMovement(nCoast);         
         WritePolygonPotentialErosion(nCoast);
      }
      
      // Now route actually-eroded sand/coarse sediment to adjacent polygons, or off-grid. sort first
      vector<vector<int> > nVVPolyAndAdjacent;
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         CGeomCoastPolygon const* pPoly = m_VCoast[nCoast].pGetPolygon(nPoly);

         vector<int> nVPolyAndAdj;

         // The first array item is the polygon coast ID
         nVPolyAndAdj.push_back(pPoly->nGetCoastID());

         if (pPoly->bDownCoastThisIter())
         {
            // Sediment is leaving this polygon in a down-coast direction: so set this as the second array item
            nVPolyAndAdj.push_back(true);

            // Set subsequent array items to be the IDs of adjacent polygons
            for (int nAdj = 0; nAdj < pPoly->nGetNumDownCoastAdjacentPolygons(); nAdj++)
            {
               int nAdjPolyID = pPoly->nGetDownCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }
         else
         {
            // Sediment is leaving this polygon in an up-coast direction: so set this as the second array item
            nVPolyAndAdj.push_back(false);

            // Set subsequent array items to be the IDs of adjacent polygons
            for (int nAdj = 0; nAdj < pPoly->nGetNumUpCoastAdjacentPolygons(); nAdj++)
            {
               int nAdjPolyID = pPoly->nGetUpCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }

         nVVPolyAndAdjacent.push_back(nVPolyAndAdj);
      }

      // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      //    WritePolygonUnsortedSequence(nCoast, nVVPolyAndAdjacent);

      // OK, now sort the array using bPolygonAndAdjCompare(), so that 'target' polygons are processed after 'source' polygons
      stable_sort(nVVPolyAndAdjacent.begin(), nVVPolyAndAdjacent.end(), bPolygonAndAdjCompare);
      
      // And check for circularities i.e. where poly X -> poly Y -> poly X. Note that we only look for two-way circularities. i.e. we ignore poly A -> poly B -> Poly C -> poly A patterns. These are probably pretty rare, however
      vector<int> VnSourcePolygons;
      for (int nPoly = 0; nPoly < static_cast<int>(nVVPolyAndAdjacent.size()); nPoly++)
      {
         for (int m = 0; m < static_cast<int>(nVVPolyAndAdjacent[nPoly].size()); m++)
         {
            if (m == 0) 
               VnSourcePolygons.push_back(nVVPolyAndAdjacent[nPoly][m]);
            
            else if (m > 1)
            {
               // Check for circularities
               vector<int>::iterator it = find(VnSourcePolygons.begin(), VnSourcePolygons.end(), nVVPolyAndAdjacent[nPoly][m]);
               if (it != VnSourcePolygons.end())
               {
                  // Uh-oh: this polygon is in the list of previously-processed source polygons. So store the numbers of the polygons with circularity
                  CGeomCoastPolygon* pPoly = m_VCoast[nCoast].pGetPolygon(nVVPolyAndAdjacent[nPoly][0]);
                  pPoly->AddCircularity(*it);
                  
                  pPoly = m_VCoast[nCoast].pGetPolygon(*it);
                  pPoly->AddCircularity(nVVPolyAndAdjacent[nPoly][0]);                     
               }
            }            
         }
      }
      
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         WritePolygonSortedSequence(nCoast, nVVPolyAndAdjacent);
      
      int nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

      // Now process all polygons and do the actual (supply-limited) unconsolidated sediment movement
      for (int n = 0; n < nNumPolygons; n++)
      {
         int nPoly = nVVPolyAndAdjacent[n][0];
         
         // Do deposition first: does this polygon have coarse deposition?
         double dCoarseDepositionTarget = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDepositionUnconsCoarse();
         
         if (dCoarseDepositionTarget > 0)
         {
            // Yes, so do deposition of coarse sediment: calculate a net increase in depth of coarse-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some coarse-sized sediment erosion) however
            double dCoarseDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, nPoly, TEXTURE_COARSE, dCoarseDepositionTarget, dCoarseDeposited);
            if (nRet != RTN_OK)
               return nRet; 
            
            double dCoarseNotDeposited = dCoarseDepositionTarget - dCoarseDeposited;
            if (dCoarseNotDeposited > 0)
            {
               m_dThisIterDepositionCoarseDiff += dCoarseNotDeposited;
            }            
         }
         
         // Does this polygon have sand deposition?
         double dSandDepositionTarget = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDepositionUnconsSand();
         
         if (dSandDepositionTarget > 0)
         {
            // Yes, so do deposition of sand sediment: calculate a net increase in depth of sand-sized unconsolidated sediment on the cells within the polygon. Note that some cells may decrease in elevation (i.e. have some sand-sized sediment erosion) however
            double dSandDeposited = 0;
            nRet = nDoUnconsDepositionOnPolygon(nCoast, nPoly, TEXTURE_SAND, dSandDepositionTarget, dSandDeposited);
            if (nRet != RTN_OK)
               return nRet; 
            
            double dSandNotDeposited = dSandDepositionTarget - dSandDeposited;
            if (dSandNotDeposited > 0)
            {
               m_dThisIterDepositionSandDiff += dSandNotDeposited;
            }            
         }
         
         // Now do erosion
         double dPotentialErosion = -m_VCoast[nCoast].pGetPolygon(nPoly)->dGetPotentialErosion();
         if (dPotentialErosion > 0)
         {
            // There is some erosion on this polygon: process this in the sequence fine, sand, coarse. Is there any fine sediment on this polygon?
            double dExistingUnconsFine = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetStoredUnconsFine();
            
            if (dExistingUnconsFine > 0)
            {
               // Yes there is, so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double dFinePotentialErosion = dPotentialErosion * m_dFineErodibilityNormalized;
                  
               // Now reduce this further, by considering the total depth of fine sediment on the polygon
               double dFineErosionTarget = tMin(dFinePotentialErosion, dExistingUnconsFine);
                              
               // OK, do the supply-limited erosion of fine sediment
               double dFineEroded = 0;
               nRet = nDoUnconsErosionOnPolygon(nCoast, nPoly, TEXTURE_FINE, dFineErosionTarget, dFineEroded);
               if (nRet != RTN_OK)
                  return nRet;
               
               if (dFineEroded > 0)
               {
                  // We eroded some fine sediment, so add to the this-iteration toal. Note that total this gets added in to the suspended load elsewhere, so no need to do it here
                  m_dThisIterBeachErosionFine += dFineEroded;
               
                  // Store the amount eroded for this polygon
                  m_VCoast[nCoast].pGetPolygon(nPoly)->SetErosionUnconsFine(-dFineEroded);
               }
            }
            
            // Is there any sand-sized sediment on this polygon?
            double dExistingUnconsSand = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetStoredUnconsSand();
            double dSandEroded = 0;            
            if (dExistingUnconsSand > 0)
            {
               // There is: so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double dSandPotentialErosion = dPotentialErosion * m_dSandErodibilityNormalized;
               
               // Now reduce this further, by considering the total depth of sand sediment on the polygon
               double dSandErosionTarget = tMin(dSandPotentialErosion, dExistingUnconsSand);
                  
               // If we have some sand sediment which we were unable to deposit on the previously-processed polygon (which could be the last-processed polygon of the previous timestep), then add this in
               if (m_dThisIterDepositionSandDiff > 0)
               {
                  dSandErosionTarget += m_dThisIterDepositionSandDiff;
                  m_dThisIterDepositionSandDiff = 0;            
               }
               
               // OK, do the supply-limited erosion of sand sediment
               nRet = nDoUnconsErosionOnPolygon(nCoast, nPoly, TEXTURE_SAND, dSandErosionTarget, dSandEroded);
               if (nRet != RTN_OK)
                  return nRet;
               
               if (dSandEroded > 0)
               {
                  // We eroded some sand sediment, so add to the this-iteration toal
                  m_dThisIterBeachErosionSand += dSandEroded;
               
                  // Store the amount eroded for this polygon
                  m_VCoast[nCoast].pGetPolygon(nPoly)->SetErosionUnconsSand(-dSandEroded);
               }               
            }
            
            // Is there any coarse sediment on this polygon?
            double dExistingUnconsCoarse = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetStoredUnconsCoarse();
            double dCoarseEroded = 0;
            if (dExistingUnconsCoarse > 0)
            {
               // There is: so crudely partition this potential value for this size class by erodibility, the result will almost always be much greater than actual (supply limited) erosion
               double dCoarsePotentialErosion = dPotentialErosion * m_dCoarseErodibilityNormalized;
               
               // Now reduce this further, by considering the total depth of coarse sediment on the polygon
               double dCoarseErosionTarget = tMin(dCoarsePotentialErosion, dExistingUnconsCoarse);
               
               // If we have some coarse sediment which we were unable to deposit on the previously-processed polygon (which could be the last-processed polygon of the previous timestep), then add this in
               if (m_dThisIterDepositionCoarseDiff > 0)
               {
                  dCoarseErosionTarget += m_dThisIterDepositionCoarseDiff;
                  m_dThisIterDepositionCoarseDiff = 0;            
               }
            
               // OK, do the supply-limited erosion of coarse sediment
               nRet = nDoUnconsErosionOnPolygon(nCoast, nPoly, TEXTURE_COARSE, dCoarseErosionTarget, dCoarseEroded);
               if (nRet != RTN_OK)
                  return nRet;

               if (dCoarseEroded > 0)
               {
                  // We eroded some coarse sediment, so add to the this-iteration toal
                  m_dThisIterBeachErosionCoarse += dCoarseEroded;
               
                  // Store the amount eroded for this polygon
                  m_VCoast[nCoast].pGetPolygon(nPoly)->SetErosionUnconsCoarse(-dCoarseEroded);
               }               
            }
            
            // OK we now have the actual values of sediment eroded from this polygon, so next determine where this eroded sand and coarse sediment goes (have to consider fine sediment too, because this goes off-grid on grid-edge polygons). Only do this if some sand or coarse was eroded on this polygon
            if ((dSandEroded + dCoarseEroded) > 0)        
            {
               CGeomCoastPolygon const* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

               if (pPolygon->bDownCoastThisIter())
               {
                  // Moving eroded sediment down-coast
                  int nNumAdjPoly = pPolygon->nGetNumDownCoastAdjacentPolygons();
                  for (int nn = 0; nn < nNumAdjPoly; nn++)
                  {
                     int nAdjPoly = pPolygon->nGetDownCoastAdjacentPolygon(nn);

                     // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nPoly << " moves sediment down-coast to polygon " << nAdjPoly << endl;

                     if (nAdjPoly == INT_NODATA)
                     {
                        // This polygon is at the grid edge
                        if (nPoly == 0)
                        {
                           // This is the polygon at the up-coast end of the coastline: uncons sediment movement is down-coast but there is no adjacent polygon!
                           if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                              LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): " << ERR << "when adjusting sediment export. Polygon " << nPoly << " is at the up-coast end of the coastline, actual sediment movement is DOWN-COAST. But there is no adjacent coast-end polygon!" << endl;
                        }
                        else if (nPoly == nNumPolygons-1)
                        {
                           // This is the polygon at the down-coast end of the coastline, and uncons sediment movement is down-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
                           if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
                           {
                              // Closed grid edges: no uncons sediment moves off-grid, nothing is removed from this polygon, so cannot adjust sediment export
                              if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                                 LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): when adjusting sediment export, polygon " << nPoly << " is at the down-coast end of the coastline, and actual sediment movement is DOWN-COAST. Since grid edges are closed, no sand or coarse unconsolidated sediment goes off-grid so cannot adjust sediment export" << endl;
                           }

                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
                           {
                              // Open grid edges, so this sediment goes off-grid
                              m_dThisIterLeftGridUnconsSand += dSandEroded;
                              m_dThisIterLeftGridUnconsCoarse += dCoarseEroded;
                           }

                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
                           {
                              // Re-circulating grid edges, so adjust the sediment exported to the polygon at the up-coast end of this coastline
                              int nOtherEndPoly = 0;
                              CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);
                              
                              if (dSandEroded > 0)
                              {
                                 pOtherEndPoly->AddDepositionUnconsSand(dSandEroded);
                              }

                              if (dCoarseEroded > 0)
                              {
                                 pOtherEndPoly->AddDepositionUnconsCoarse(dCoarseEroded);
                              }
                           }
                        }
                     }
                     else
                     {
                        // This polygon is not at the grid edge
                        CGeomCoastPolygon* pAdjPolygon = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
                        double dBoundaryShare = pPolygon->dGetDownCoastAdjacentPolygonBoundaryShare(nn);

                        if (dSandEroded > 0)
                        {
                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nAdjPoly << ", dDeltaSand was = " << pAdjPolygon->dGetDepositionUnconsSand() << " m_dThisIterErosionSandDiff = " << m_dThisIterErosionSandDiff;

                           pAdjPolygon->AddDepositionUnconsSand(dSandEroded * dBoundaryShare);

                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << " dDeltaSand now = " << pAdjPolygon->dGetDepositionUnconsSand() << endl;
                        }

                        if (dCoarseEroded > 0)
                        {
                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nAdjPoly << ", dDeltaCoarse was = " << pAdjPolygon->dGetDepositionUnconsCoarse() << " m_dThisIterErosionCoarseDiff = " << m_dThisIterErosionCoarseDiff;

                           pAdjPolygon->AddDepositionUnconsCoarse(dCoarseEroded * dBoundaryShare);

                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << " dDeltaCoarse now = " << pAdjPolygon->dGetDepositionUnconsCoarse() << endl;
                        }
                     }
                  }

                  // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  //    LogStream << " dFineEroded = " << dFineEroded << " dSandEroded = " << dSandEroded << " dCoarseEroded = " << dCoarseEroded << endl;
               }
               else
               {
                  // Moving eroded sediment up-coast
                  int nNumAdjPoly = pPolygon->nGetNumUpCoastAdjacentPolygons();
                  for (int nn = 0; nn < nNumAdjPoly; nn++)
                  {
                     int nAdjPoly = pPolygon->nGetUpCoastAdjacentPolygon(nn);
                     // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nPoly << " moves sediment up-coast to polygon " << nAdjPoly << endl;

                     if (nAdjPoly == INT_NODATA)
                     {
                        // This polygon is at the grid edge
                        if (nPoly == nNumPolygons-1)
                        {
                           // This is the polygon at the down-coast end of the coastline: uncons sediment movement is up-coast but there is no adjacent polygon!
                           if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                              LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): " << ERR << "when adjusting sediment export. Polygon " << nPoly << " is at the down-coast end of the coastline, actual sediment movement is UP-COAST. But there is no adjacent coast-end polygon!" << endl;
                        }
                        else if (nPoly == 0)
                        {
                           // This is the polygon at the up-coast end of the coastline, and uncons sediment movement is up-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
                           if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
                           {
                              // Closed grid edges: no uncons sediment moves off-grid, nothing is removed from this polygon, so cannot adjust sediment export
                              if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                                 LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): when adjusting sediment export, polygon " << nPoly << " is at the up-coast end of the coastline, and actual sediment movement is UP-COAST. Since grid edges are closed, no sand or coarse unconsolidated sediment goes off-grid so cannot adjust sediment export" << endl;
                           }

                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
                           {
                              // Open grid edges, so this sediment goes off-grid
                              m_dThisIterLeftGridUnconsSand += dSandEroded;
                              m_dThisIterLeftGridUnconsCoarse += dCoarseEroded;
                           }

                           else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
                           {
                              // Re-circulating grid edges, so adjust the sediment exported to the polygon at the up-coast end of this coastline TODO 016 Check whether this causes mass balance problems, depending on the sequence of polygon processing
                              int nOtherEndPoly = 0;
                              CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);

                              if (dSandEroded > 0)
                              {
                                 pOtherEndPoly->AddDepositionUnconsSand(dSandEroded);
                              }

                              if (dCoarseEroded > 0)
                              {
                                 pOtherEndPoly->AddDepositionUnconsCoarse(dCoarseEroded);
                              }
                           }
                        }
                     }
                     else
                     {
                        // This polygon is not at the grid edge
                        CGeomCoastPolygon* pAdjPolygon = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
                        double dBoundaryShare = pPolygon->dGetUpCoastAdjacentPolygonBoundaryShare(nn);

                        if (dSandEroded > 0)
                        {
                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nAdjPoly << ", dDeltaSand was = " << pAdjPolygon->dGetDepositionUnconsSand() << " m_dThisIterErosionSandDiff = " << m_dThisIterErosionSandDiff;

                           pAdjPolygon->AddDepositionUnconsSand(dSandEroded * dBoundaryShare);

                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << " dDeltaSand now = " << pAdjPolygon->dGetDepositionUnconsSand() << endl;
                        }

                        if (dCoarseEroded > 0)
                        {
                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << "Timestep " << m_ulIter << " (" << strDispSimTime(m_dSimElapsed) << "): polygon " << nAdjPoly << ", dDeltaCoarse was = " << pAdjPolygon->dGetDepositionUnconsCoarse() << " m_dThisIterErosionCoarseDiff = " << m_dThisIterErosionCoarseDiff;

                           pAdjPolygon->AddDepositionUnconsCoarse(+dCoarseEroded * dBoundaryShare);

                           // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                           //    LogStream << " dDeltaCoarse now = " << pAdjPolygon->dGetDepositionUnconsCoarse() << endl;
                        }
                     }
                  }
               }      
            }
            
            // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            //    LogStream << " dFineEroded = " << dFineEroded << " dSandEroded = " << dSandEroded << " dCoarseEroded = " << dCoarseEroded << endl;
            
         }     // if (dPotentialErosion > 0)
      }     // for (int n = 0; n < nNumPolygons; n++)
      
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         WritePolygonActualMovement(nCoast, nVVPolyAndAdjacent);      
   }
   
   return RTN_OK;
}

//===============================================================================================================================
//! Update the values of pre-existing unconsolidated sediment, for all three size classes, to include unconsolidated sediment derived from platform erosion and/or cliff collapse
//===============================================================================================================================
void CSimulation::AllPolygonsUpdateStoredUncons(int const nCoast)
{
   int nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

   // Update the polygons, unconsolidated sand and coarse only (any fine sediment from platform erosion and cliff collapse goes to suspension)
   for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
   {
      CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
      
      double dSand = pPolygon->dGetStoredUnconsSand() + pPolygon->dGetUnconsSandFromShorePlatform() + pPolygon->dGetCliffCollapseUnconsSandDeposition();      
      pPolygon->SetStoredUnconsSand(dSand);
      
      double dCoarse = pPolygon->dGetStoredUnconsCoarse() + pPolygon->dGetUnconsCoarseFromShorePlatform() + pPolygon->dGetCliffCollapseUnconsCoarseDeposition();
      pPolygon->SetStoredUnconsCoarse(dCoarse);
   }   
}
