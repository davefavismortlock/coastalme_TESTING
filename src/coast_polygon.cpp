/*!
 *
 * \file coast_polygon.cpp
 * \brief CGeomCoastPolygon routines
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <assert.h>
#include <iostream>
using std::cerr;

#include "cme.h"
#include "coast_polygon.h"

//! Constructor with 8 parameters and initialization list
CGeomCoastPolygon::CGeomCoastPolygon(int const nGlobalID, int const nCoastID, int const nNode, int const nProfileUpCoast, int const nProfileDownCoast, vector<CGeom2DPoint> const* pVIn, int const nLastPointUpCoast, const int nLastPointDownCoast, CGeom2DIPoint const* PtiNode, CGeom2DIPoint const* PtiAntinode, int const nPointInPolygonStartPoint)
:
//    m_bIsPointedSeaward(true),
   m_bDownCoastThisIter(false),
   m_nGlobalID(nGlobalID),
   m_nCoastID(nCoastID),
   m_nCoastNode(nNode),
   m_nProfileUpCoast(nProfileUpCoast),
   m_nProfileDownCoast(nProfileDownCoast),
   m_nProfileUpCoastNumPointsUsed(nLastPointUpCoast),
   m_nProfileDownCoastNumPointsUsed(nLastPointDownCoast),
   m_nNumCells(0),
   m_nPointInPolygonSearchStartPoint(nPointInPolygonStartPoint),   
   m_dAvgUnconsD50(0),   
   m_dSeawaterVolume(0),
   m_dPotentialErosionAllUncons(0),
   m_dErosionUnconsFine(0),
   m_dErosionUnconsSand(0),
   m_dErosionUnconsCoarse(0),
   m_dDepositionUnconsFine(0),
   m_dDepositionUnconsSand(0),
   m_dDepositionUnconsCoarse(0),
   m_dCliffCollapseErosionFine(0),
   m_dCliffCollapseErosionSand(0),
   m_dCliffCollapseErosionCoarse(0),
   m_dCliffCollapseTalusSand(0),
   m_dCliffCollapseTalusCoarse(0),
   m_dSandFromPlatformErosion(0),
   m_dCoarseFromPlatformErosion(0),
   m_dStoredUnconsFine(0),
   m_dStoredUnconsSand(0),
   m_dStoredUnconsCoarse(0),
   m_PtiNode(*PtiNode),
   m_PtiAntinode(*PtiAntinode)
{
   m_VPoints = *pVIn;
}

//! Destructor
CGeomCoastPolygon::~CGeomCoastPolygon(void)
{
}

// void CGeomCoastPolygon::SetNotPointed(void)
// {
//    m_bIsPointedSeaward = false;
// }
//
// bool CGeomCoastPolygon::bIsPointed(void) const
// {
//    return m_bIsPointedSeaward;
// }

//! Set a flag to say whether sediment movement on this polygon is downcoast this iteration
void CGeomCoastPolygon::SetDownCoastThisIter(bool const bFlag)
{
   m_bDownCoastThisIter = bFlag;
}

//! Is sediment movement on this polygon downcoast this iteration?
bool CGeomCoastPolygon::bDownCoastThisIter(void) const
{
   return m_bDownCoastThisIter;
}

//! Get the global ID
int CGeomCoastPolygon::nGetGlobalID(void) const
{
   return m_nGlobalID;
}

//! Get the coast ID
int CGeomCoastPolygon::nGetCoastID(void) const
{
   return m_nCoastID;
}

// void CGeomCoastPolygon::SetCoastNode(int const nNode)
// {
//    m_nCoastNode = nNode;
// }

//! Get the coast node point
int CGeomCoastPolygon::nGetNodeCoastPoint(void) const
{
   return m_nCoastNode;
}

//! Get the grid co-ordinates of the cell on which the node sits
CGeom2DIPoint* CGeomCoastPolygon::pPtiGetNode(void)
{
   return &m_PtiNode;

}

//! Get the anti-node (raster-grid CRS) which is at other (seaward) end of the polygon from the node
CGeom2DIPoint* CGeomCoastPolygon::pPtiGetAntiNode(void)
{
   return &m_PtiAntinode;
}

//! Set the number of cells in the polygon
void CGeomCoastPolygon::SetNumCellsInPolygon(int const nCells)
{
   m_nNumCells = nCells;
}

// //! Get the number of cells in the polygon
// int CGeomCoastPolygon::nGetNumCellsinPolygon(void) const
// {
//    return m_nNumCells;
// }

//! Return the number of the up-coast profile
int CGeomCoastPolygon::nGetUpCoastProfile(void) const
{
   return m_nProfileUpCoast;
}

//! Return the number of the down-coast profile
int CGeomCoastPolygon::nGetDownCoastProfile(void) const
{
   return m_nProfileDownCoast;
}

// void CGeomCoastPolygon::SetBoundary(vector<CGeom2DPoint> const* pVIn)
// {
//    m_VPoints = *pVIn;
// }

// vector<CGeom2DPoint>* CGeomCoastPolygon::pPtVGetBoundary(void)
// {
//    return &m_VPoints;
// }

//! Get the co-ordinates (external CRS) of a specified point on the polygon's boundary
CGeom2DPoint* CGeomCoastPolygon::pPtGetBoundaryPoint(int const nPoint)
{
   // TODO 055 No check to see if nPoint < m_VPoints.size()
   return &m_VPoints[nPoint];
}

//! Get the number of points in the polygon's boundary
int CGeomCoastPolygon::nGetBoundarySize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

//! Return the number of points in the up-coast profile
int CGeomCoastPolygon::nGetUpCoastProfileNumPointsUsed(void) const
{
   return m_nProfileUpCoastNumPointsUsed;
}

//! Return the number of points in the down-coast profile
int CGeomCoastPolygon::nGetDownCoastProfileNumPointsUsed(void) const
{
   return m_nProfileDownCoastNumPointsUsed;
}

//! Set the volume of seawater in the coast polygon
void CGeomCoastPolygon::SetSeawaterVolume(const double dDepth)
{
   m_dSeawaterVolume = dDepth;
}

//! Get the volume of seawater in the coast polygon
double CGeomCoastPolygon::dGetSeawaterVolume(void) const
{
   return m_dSeawaterVolume;
}

//! Adds in potential erosion of unconsolidated sediment (all size classes) on this polygon (m_dPotentialErosionAllUncons is <= 0)
void CGeomCoastPolygon::AddPotentialErosion(double const dDepth)
{
   m_dPotentialErosionAllUncons += dDepth;
}

//! Returns this timestep's total change in depth of unconsolidated sediment (all size classes) on this polygon (m_dPotentialErosionAllUncons <= 0)
double CGeomCoastPolygon::dGetPotentialErosion(void) const
{
   return m_dPotentialErosionAllUncons;
}

//! Sets a value (must be < 0) for this timestep's erosion of fine unconsolidated sediment on this polygon
void CGeomCoastPolygon::SetErosionUnconsFine(double const dDepth)
{
   m_dErosionUnconsFine = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of fine unconsolidated sediment on this polygon
double CGeomCoastPolygon::dGetErosionUnconsFine(void) const
{
   return m_dErosionUnconsFine;
}

//! Sets a value (must be < 0) for this timestep's erosion of sand-sized unconsolidated sediment on this polygon
void CGeomCoastPolygon::SetErosionUnconsSand(double const dDepth)
{
   m_dErosionUnconsSand = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of sand-sized unconsolidated sediment on this polygon
double CGeomCoastPolygon::dGetErosionUnconsSand(void) const
{
   return m_dErosionUnconsSand;
}

//! Sets a value (must be < 0) for this timestep's erosion of coarse unconsolidated sediment on this polygon
void CGeomCoastPolygon::SetErosionUnconsCoarse(double const dDepth)
{
   m_dErosionUnconsCoarse = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of coarse unconsolidated sediment on this polygon
double CGeomCoastPolygon::dGetErosionUnconsCoarse(void) const
{
   return m_dErosionUnconsCoarse;
}

//! Returns this timestep's total (all size classes) eroson of unconsolidated sediment on this polygon, as a -ve depth in m
double CGeomCoastPolygon::dGetErosionAllUncons(void) const
{
   return m_dErosionUnconsFine + m_dErosionUnconsSand + m_dErosionUnconsCoarse;
}

// //! Re-initializes this timestep's eposition of unconsolidated fine sediment on this polygon
// void CGeomCoastPolygon::SetZeroDepositionUnconsFine(void)
// {
//    m_dDepositionUnconsFine = 0;
// }

//! Returns this timestep's deposition of fine unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetDepositionUnconsFine(void) const
{
   return m_dDepositionUnconsFine;
}

// //! Re-initializes this timestep's deposition of unconsolidated sand sediment on this polygon
// void CGeomCoastPolygon::SetZeroDepositionUnconsSand(void)
// {
//    m_dDepositionUnconsSand = 0;
// }

//! Adds a depth (in m) of sand-sized unconsolidated sediment to this timestep's deposition of unconsolidated coarse sediment on this polygon
void CGeomCoastPolygon::AddDepositionUnconsSand(double const dDepth)
{
   m_dDepositionUnconsSand += dDepth;
}

//! Returns this timestep's deposition of sand-sized unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetDepositionUnconsSand(void) const
{
   return m_dDepositionUnconsSand;
}

// //! Re-initializes this timestep's deposition of unconsolidated coarse sediment on this polygon
// void CGeomCoastPolygon::SetZeroDepositionUnconsCoarse(void)
// {
//    m_dDepositionUnconsCoarse = 0;
// }

//! Adds a depth (in m) of coarse unconsolidated sediment to this timestep's deposition of unconsolidated coarse sediment on this polygon (+ve)
void CGeomCoastPolygon::AddDepositionUnconsCoarse(double const dDepth)
{
   m_dDepositionUnconsCoarse += dDepth;
}

//! Returns this timestep's deposition of coarse unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetDepositionUnconsCoarse(void) const
{
   return m_dDepositionUnconsCoarse;
}

//! Returns this timestep's total (all size classes) deposition of unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetDepositionAllUncons(void) const
{
   return m_dDepositionUnconsFine + m_dDepositionUnconsSand + m_dDepositionUnconsCoarse;
}

//! Sets all up-coast adjacent polygons
void CGeomCoastPolygon::SetUpCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnUpCoastAdjacentPolygon = *pnVPolygons;
}

//! Gets a single up-coast adjacent polygon
int CGeomCoastPolygon::nGetUpCoastAdjacentPolygon(int const nIndex) const
{
//    assert(nIndex < m_VnUpCoastAdjacentPolygon.size());
   return m_VnUpCoastAdjacentPolygon[nIndex];
}

//! Gets all up-coast adjacent polygons
int CGeomCoastPolygon::nGetNumUpCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnUpCoastAdjacentPolygon.size());
}

//! Sets all down-coast adjacent polygons
void CGeomCoastPolygon::SetDownCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnDownCoastAdjacentPolygon = *pnVPolygons;
}

//! Gets a single down-coast adjacent polygon
int CGeomCoastPolygon::nGetDownCoastAdjacentPolygon(int const nIndex) const
{
//    assert(nIndex < m_VnDownCoastAdjacentPolygon.size());
   return m_VnDownCoastAdjacentPolygon[nIndex];
}

//! Gets all down-coast adjacent polygons
int CGeomCoastPolygon::nGetNumDownCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnDownCoastAdjacentPolygon.size());
}

//! Sets the boundary shares for all up-coast adjacent polygons
void CGeomCoastPolygon::SetUpCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdUpCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

//! Gets the boundary shares for all up-coast adjacent polygons
double CGeomCoastPolygon::dGetUpCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // TODO 055 No check to see if nIndex < m_VdUpCoastAdjacentPolygonBoundaryShare.size()
   return m_VdUpCoastAdjacentPolygonBoundaryShare[nIndex];
}

//! Sets the boundary shares for all down-coast adjacent polygons
void CGeomCoastPolygon::SetDownCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdDownCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

//! Gets the boundary shares for all down-coast adjacent polygons
double CGeomCoastPolygon::dGetDownCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // TODO 055 No check to see if nIndex < m_VdDownCoastAdjacentPolygonBoundaryShare.size()
   return m_VdDownCoastAdjacentPolygonBoundaryShare[nIndex];
}

//! Returns the start point for a point-in-polygon search
int CGeomCoastPolygon::nGetPointInPolygonSearchStartPoint(void) const
{
   return m_nPointInPolygonSearchStartPoint;
}

//! Set the average d50 for unconsolidated sediment in this polygon
void CGeomCoastPolygon::SetAvgUnconsD50(double const dD50)
{
   m_dAvgUnconsD50 = dD50;
}

//! Get the average d50 for unconsolidated sediment in this polygon
double CGeomCoastPolygon::dGetAvgUnconsD50(void) const
{
   return m_dAvgUnconsD50;
}

//! Instantiates the pure virtual function in the abstract parent class, so that CGeomCoastPolygon is not an abstract class
void CGeomCoastPolygon::Display(void)
{
}

//! Add a circularity to this polygon
void CGeomCoastPolygon::AddCircularity(int const nPoly)
{
   m_VnCircularityWith.push_back(nPoly);
}

//! Get all circularities for this polygon
vector<int> CGeomCoastPolygon::VnGetCircularities(void) const
{
   return m_VnCircularityWith;
}

//! Add to the this-iteration total of unconsolidated fine sediment from cliff collapse in this polygon
void CGeomCoastPolygon::AddCliffCollapseErosionFine(double const dDepth)
{
   m_dCliffCollapseErosionFine += dDepth;
}

//! Get the this-iteration total of unconsolidated fine sediment from cliff collapse in this polygon
double CGeomCoastPolygon::dGetCliffCollapseErosionFine(void) const
{
   return m_dCliffCollapseErosionFine;
}

//! Add to the this-iteration total of unconsolidated sand sediment from cliff collapse in this polygon
void CGeomCoastPolygon::AddCliffCollapseErosionSand(double const dDepth)
{
   m_dCliffCollapseErosionSand += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment from cliff collapse in this polygon
double CGeomCoastPolygon::dGetCliffCollapseErosionSand(void) const
{
   return m_dCliffCollapseErosionSand;
}

//! Add to the this-iteration total of unconsolidated coarse sediment from cliff collapse in this polygon
void CGeomCoastPolygon::AddCliffCollapseErosionCoarse(double const dDepth)
{
   m_dCliffCollapseErosionCoarse += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment from cliff collapse in this polygon
double CGeomCoastPolygon::dGetCliffCollapseErosionCoarse(void) const
{
   return m_dCliffCollapseErosionCoarse;
}

//! Add to the this-iteration total of unconsolidated sand sediment deposited from cliff collapse in this polygon
void CGeomCoastPolygon::AddCliffCollapseUnconsSandDeposition(double const dDepth)
{
   m_dCliffCollapseTalusSand += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment deposited from cliff collapse in this polygon
double CGeomCoastPolygon::dGetCliffCollapseUnconsSandDeposition(void) const
{
   return m_dCliffCollapseTalusSand;
}

//! Add to the this-iteration total of unconsolidated coarse sediment deposited from cliff collapse in this polygon
void CGeomCoastPolygon::AddCliffCollapseUnconsCoarseDeposition(double const dDepth)
{  
   m_dCliffCollapseTalusCoarse += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment deposited from cliff collapse in this polygon
double CGeomCoastPolygon::dGetCliffCollapseUnconsCoarseDeposition(void) const
{
   return m_dCliffCollapseTalusCoarse;
}

//! Add to the this-iteration total of unconsolidated sand sediment derived from shore platform erosion in this polygon
void CGeomCoastPolygon::AddUnconsSandFromShorePlatform(double const dDepth)
{
   m_dSandFromPlatformErosion += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment derived from shore platform erosion in this polygon
double CGeomCoastPolygon::dGetUnconsSandFromShorePlatform(void) const
{
   return m_dSandFromPlatformErosion;
}

//! Add to the this-iteration total of unconsolidated coarse sediment derived from shore platform erosion in this polygon
void CGeomCoastPolygon::AddUnconsCoarseFromShorePlatform(double const dDepth)
{
   m_dCoarseFromPlatformErosion += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment derived from shore platform erosion in this polygon
double CGeomCoastPolygon::dGetUnconsCoarseFromShorePlatform(void) const
{
   return m_dCoarseFromPlatformErosion;
}

//! Set the value of stored unconsolidated fine sediment stored in this polygon
void CGeomCoastPolygon::SetStoredUnconsFine(double const dDepth)
{
   m_dStoredUnconsFine = dDepth;
}

//! Get the value of stored unconsolidated fine sediment stored in this polygon
double CGeomCoastPolygon::dGetStoredUnconsFine(void) const
{
   return m_dStoredUnconsFine;
}

//! Set the value of stored unconsolidated sand sediment stored in this polygon
void CGeomCoastPolygon::SetStoredUnconsSand(double const dDepth)
{
   m_dStoredUnconsSand = dDepth;
}

//! Get the value of stored unconsolidated sand sediment stored in this polygon
double CGeomCoastPolygon::dGetStoredUnconsSand(void) const
{
   return m_dStoredUnconsSand;
}

//! Set the value of stored unconsolidated coarse sediment stored in this polygon
void CGeomCoastPolygon::SetStoredUnconsCoarse(double const dDepth)
{
   m_dStoredUnconsCoarse = dDepth;
}

//! Get the value of stored unconsolidated coarse sediment stored in this polygon
double CGeomCoastPolygon::dGetStoredUnconsCoarse(void) const
{
   return m_dStoredUnconsCoarse;
}

