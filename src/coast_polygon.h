/*!
 * \class CGeomCoastPolygon
 * \brief Geometry class used for coast polygon objects
 * \details TODO 001 This is a more detailed description of the CRWCoast class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
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
   //! Is the movement of unconsolidated sediment on this polygon down-coast during this iteration?
   bool m_bDownCoastThisIter;

   // Does the polygon meet at a point at its seaward end? (is it roughly triangular?)
//   bool m_bIsPointedSeaward,

   //! The simulation-global number of this polygon
   int m_nGlobalID;

   //! This-coast-only number of this polygon
   int m_nCoastID;

   //! The point on this polygon's coastline segment with maximum concave curvature, roughly at the middle of the coastline segment
   int m_nCoastNode;

   //! The normal profile which bounds the polygon in the up-coast direction
   int m_nProfileUpCoast;

   //! The normal profile which bounds the polygon in the down-coast direction
   int m_nProfileDownCoast;

   //! The number of points from the up-coast normal which are part of this polygon (less than the normal's full length if the polygon is triangular)
   int m_nProfileUpCoastNumPointsUsed;

   //! The number of points from the down-coast normal which are part of this polygon (less than the normal's full length if the polygon is triangular)
   int m_nProfileDownCoastNumPointsUsed;

   //! The number of cells in the polygon
   int m_nNumCells;

   //! The number of the vector point from which we start the point-in-polygon search
   int m_nPointInPolygonSearchStartPoint;
      
   // Note: all sediment depths are in m, and here cover the area of a single raster cell: to convert to a volume, multiply by m_dCellArea

   //! The average d50 of unconsolidated sediment on this polygon
   double m_dAvgUnconsD50;

   //! The volume (m3) of seawater within the polygon
   double m_dSeawaterVolume;

   //! Potential (ignoring supply-limitation) erosion (all size classes) as a depth during this timestep (-ve)
   double m_dPotentialErosionAllUncons;

   //! Erosion (considering supply-limitation) of fine-sized sediment as a depth this timestep (-ve)
   double m_dErosionUnconsFine;

   //! Erosion (considering supply-limitation) of sand-sized sediment as a depth this timestep (-ve)
   double m_dErosionUnconsSand;

   //! Erosion (considering supply-limitation) of coarse-sized sediment as a depth this timestep (-ve)
   double m_dErosionUnconsCoarse;

   //! Deposition of fine-sized sediment as a depth this timestep (+ve)
   double m_dDepositionUnconsFine;

   //! Deposition of sand-sized sediment as a depth this timestep (+ve)
   double m_dDepositionUnconsSand;

   //! Deposition of coarse-sized sediment as a depth this timestep (+ve)
   double m_dDepositionUnconsCoarse;

   //! Depth of eroded fine sediment from cliff collapse
   double m_dCliffCollapseErosionFine;

   //! Depth of eroded sand sediment from cliff collapse
   double m_dCliffCollapseErosionSand;

   //! Depth of eroded coarse sediment from cliff collapse
   double m_dCliffCollapseErosionCoarse;

   //! Depth of unconsolidated sand talus from cliff collapse
   double m_dCliffCollapseTalusSand;

   //! Depth of unconsolidated coarse talus from cliff collapse
   double m_dCliffCollapseTalusCoarse;

   //! Depth of unconsolidated sand sediment from shore platform erosion
   double m_dSandFromPlatformErosion;

   //! Depth of unconsolidated coarse sediment from shore platform erosion
   double m_dCoarseFromPlatformErosion;

   //! Depth of pre-existing unconsolidated fine sediment
   double m_dStoredUnconsFine;

   //! Depth of pre-existing unconsolidated sand sediment
   double m_dStoredUnconsSand;

   //! Depth of pre-existing unconsolidated coarse sedimet
   double m_dStoredUnconsCoarse;

   //! Co-ordinates of the coast node cell (raster-grid CRS)
   CGeom2DIPoint m_PtiNode;

   //! Co-ordinates of the cell (raster-grid CRS) which is at other (seaward) end of the polygon
   CGeom2DIPoint m_PtiAntinode;

   //! The ID(s) of the up-coast adjacent polygon(s)
   vector<int> m_VnUpCoastAdjacentPolygon;

   //! The ID(s) of the down-coast adjacent polygon(s)
   vector<int> m_VnDownCoastAdjacentPolygon;

   //! If this polygon has a circular unconsolidated-sediment-movement relationship with one or more other polygons, the cost-only numbers of these polygons
   vector<int> m_VnCircularityWith;

   //! The boundary share(s) (0 to 1) with adjacent up-coast polygon(s)
   vector<double> m_VdUpCoastAdjacentPolygonBoundaryShare;

   //! The boundary share(s) (0 to 1) with adjacent up-coast polygon(s)
   vector<double> m_VdDownCoastAdjacentPolygonBoundaryShare;

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
   // int nGetNumCellsinPolygon(void) const;

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

   // void SetZeroDepositionUnconsFine(void);
   double dGetDepositionUnconsFine(void) const;
   // void SetZeroDepositionUnconsSand(void);
   void AddDepositionUnconsSand(double const);
   double dGetDepositionUnconsSand(void) const;
   // void SetZeroDepositionUnconsCoarse(void);
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
   vector<int> VnGetCircularities(void) const;
   
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

