/*!
 * \class CGeomCoastPolygon
 * \brief Geometry class used for coast polygon objects
 * \details TODO This is a more detailed description of the CRWCoast class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file coast_polygon.h
 * \brief Contains CGeomCoastPolygon definitions
 *
 */

#ifndef COASTPOLYGON_H
#define COASTPOLYGON_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "2d_shape.h"

class CGeomCoastPolygon : public CA2DShape
{
private:
   bool
//       m_bIsPointedSeaward,                // Does the polygon meet at a point at its seaward end? (is it roughly triangular?)
      m_bDownCoastThisIter;

   int
      m_nGlobalID,                        // The simulation-global number of this polygon
      m_nCoastID,                         // This-coast-only number of this polygon
      m_nCoastNode,                       // The point on this polygon's coastline segment with maximum concave curvature, roughly at the middle of the coastline segment
      m_nProfileUpCoast,                  // The normal profile which bounds the polygon in the up-coast direction
      m_nProfileDownCoast,                // Ditto for the down-coast direction
      m_nProfileUpCoastNumPointsUsed,     // The number of points from the up-coast normal which are part of this polygon (less than the normal's full length if the polygon is triangular)
      m_nProfileDownCoastNumPointsUsed,   // Ditto for the down-coast normal
      m_nNumCells,                        // The number of cells in the polygon
      m_nPointInPolygonSearchStartPoint;  // The number of the vector point from which we start the point-in-polygon search
      
   // Note: all sediment depths are in m, and here cover the area of a single raster cell: to convert to a volume, multiply by m_dCellArea
   double
      m_dAvgUnconsD50,                    // The average d50 of unconsolidated sediment on this polygon
      m_dSeawaterVolume,                  // The volume (m3) of seawater within the polygon
      m_dPotentialErosionAllUncons,       // Potential (ignoring supply-limitation) erosion (all size classes) as a depth during this timestep (-ve)
      m_dErosionUnconsFine,               // Erosion (considering supply-limitation) of fine-sized sediment as a depth this timestep (-ve)
      m_dErosionUnconsSand,               // Erosion (considering supply-limitation) of sand-sized sediment as a depth this timestep (-ve)
      m_dErosionUnconsCoarse,             // Erosion (considering supply-limitation) of coarse-sized sediment as a depth this timestep (-ve)
      m_dDepositionUnconsFine,            // Deposition of fine-sized sediment as a depth this timestep (+ve)
      m_dDepositionUnconsSand,            // Deposition of sand-sized sediment as a depth this timestep (+ve)
      m_dDepositionUnconsCoarse,          // Deposition of coarse-sized sediment as a depth this timestep (+ve)
      m_dCliffCollapseErosionFine,        // Depth of eroded fine sediment from cliff collapse
      m_dCliffCollapseErosionSand,        // Depth of eroded sand sediment from cliff collapse
      m_dCliffCollapseErosionCoarse,      // Depth of eroded coarse sediment from cliff collapse
      m_dCliffCollapseTalusSand,          // Depth of unconsolidated sand talus from cliff collapse
      m_dCliffCollapseTalusCoarse,        // Depth of unconsolidated coarse talus from cliff collapse
      m_dSandFromPlatformErosion,         // Depth of unconsolidated sand sediment from shore platform erosion
      m_dCoarseFromPlatformErosion,       // Depth of unconsolidated coarse sediment from shore platform erosion
      m_dStoredUnconsFine,           // Depth of pre-existing unconsolidated fine sediment
      m_dStoredUnconsSand,           // Depth of pre-existing unconsolidated sand sediment
      m_dStoredUnconsCoarse;         // Depth of pre-existing unconsolidated coarse sediment

   CGeom2DIPoint
      m_PtiNode,                          // Co-ords of the coast node cell (raster-grid CRS)
      m_PtiAntinode;                      // Co-ords of the cell (raster-grid CRS) which is at other (seaward) end of the polygon

   vector<int>
      m_VnUpCoastAdjacentPolygon,
      m_VnDownCoastAdjacentPolygon,
      m_VnCircularityWith;                // If this polygon has a circular unconsolidated-sediment-movement relationship with one or more other polygons, the cost-only numbers of these polygons

   vector<double>
      m_VdUpCoastAdjacentPolygonBoundaryShare,
      m_VdDownCoastAdjacentPolygonBoundaryShare;

public:
   CGeomCoastPolygon(int const, int const, int const, int const, int const, vector<CGeom2DPoint> const*, int const, int const, CGeom2DIPoint const*, CGeom2DIPoint const*, int const);
   ~CGeomCoastPolygon(void);

   void SetDownCoastThisIter(bool const);
   bool bDownCoastThisIter(void) const;

   int nGetGlobalID(void) const;

   int nGetCoastID(void) const;

//    void SetCoastNode(int const);
   int nGetNodeCoastPoint(void) const;
   CGeom2DIPoint* pPtiGetNode(void);
   CGeom2DIPoint* pPtiGetAntiNode(void);

//    void SetNotPointed(void);
//    bool bIsPointed(void) const;

   void SetNumCellsInPolygon(int const);
   int nGetNumCellsinPolygon(void) const;

   int nGetUpCoastProfile(void) const;
   int nGetDownCoastProfile(void) const;

//    void SetBoundary(vector<CGeom2DPoint> const*);
//    vector<CGeom2DPoint>* pPtVGetBoundary(void);
   CGeom2DPoint* pPtGetBoundaryPoint(int const);
   int nGetBoundarySize(void) const;

   int nGetUpCoastProfileNumPointsUsed(void) const;
   int nGetDownCoastProfileNumPointsUsed(void) const;

   void SetSeawaterVolume(const double);
   double dGetSeawaterVolume(void) const;

   void AddPotentialErosion(double const);
   double dGetPotentialErosion(void) const;

   void SetErosionUnconsFine(double const);
   double dGetErosionUnconsFine(void) const;
   void SetErosionUnconsSand(double const);
   double dGetErosionUnconsSand(void) const;
   void SetErosionUnconsCoarse(double const);
   double dGetErosionUnconsCoarse(void) const;
   double dGetErosionAllUncons(void) const;

   void SetZeroDepositionUnconsFine(void);
   double dGetDepositionUnconsFine(void) const;
   void SetZeroDepositionUnconsSand(void);
   void AddDepositionUnconsSand(double const);
   double dGetDepositionUnconsSand(void) const;
   void SetZeroDepositionUnconsCoarse(void);
   void AddDepositionUnconsCoarse(double const);
   double dGetDepositionUnconsCoarse(void) const;
   double dGetDepositionAllUncons(void) const;

   void SetUpCoastAdjacentPolygons(vector<int> const*);
   int nGetUpCoastAdjacentPolygon(int const) const;
   int nGetNumUpCoastAdjacentPolygons(void) const;

   void SetDownCoastAdjacentPolygons(vector<int> const*);
   int nGetDownCoastAdjacentPolygon(int const) const;
   int nGetNumDownCoastAdjacentPolygons(void) const;

   void SetUpCoastAdjacentPolygonBoundaryShares(vector<double> const*);
   double dGetUpCoastAdjacentPolygonBoundaryShare(int const) const;

   void SetDownCoastAdjacentPolygonBoundaryShares(vector<double> const*);
   double dGetDownCoastAdjacentPolygonBoundaryShare(int const) const;

   int nGetPointInPolygonSearchStartPoint(void) const;

   void SetAvgUnconsD50(double const);
   double dGetAvgUnconsD50(void) const;

   void Display(void) override;
   
   void AddCircularity(int const);
   vector<int> VnGetCircularities(void);
   
   void AddCliffCollapseErosionFine(double const);
   double dGetCliffCollapseErosionFine(void) const;
   void AddCliffCollapseErosionSand(double const);
   double dGetCliffCollapseErosionSand(void) const;
   void AddCliffCollapseErosionCoarse(double const);
   double dGetCliffCollapseErosionCoarse(void) const;   
   
   void AddCliffCollapseUnconsSandDeposition(double const);
   double dGetCliffCollapseUnconsSandDeposition(void) const;
   void AddCliffCollapseUnconsCoarseDeposition(double const);
   double dGetCliffCollapseUnconsCoarseDeposition(void) const;   
   
   void AddUnconsSandFromShorePlatform(double const);
   double dGetUnconsSandFromShorePlatform(void) const;
   void AddUnconsCoarseFromShorePlatform(double const);
   double dGetUnconsCoarseFromShorePlatform(void) const;   
   
   void SetStoredUnconsFine(double const);
   double dGetStoredUnconsFine(void) const;
   void SetStoredUnconsSand(double const);
   double dGetStoredUnconsSand(void) const;
   void SetStoredUnconsCoarse(double const);
   double dGetStoredUnconsCoarse(void) const;
};
#endif //COASTPOLYGON_H

