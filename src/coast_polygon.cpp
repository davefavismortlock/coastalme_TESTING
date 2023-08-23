/*!
 *
 * \file coast_polygon.cpp
 * \brief CGeomCoastPolygon routines
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
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

#include "cme.h"
#include "coast_polygon.h"


//! Constructor with 8 parameters
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
   m_dStoredUnconsSand(0),           // Depth of pre-existing unconsolidated sand sediment
   m_dStoredUnconsCoarse(0),
   m_PtiNode(*PtiNode),
   m_PtiAntinode(*PtiAntinode)
{
   m_VPoints = *pVIn;
}

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


void CGeomCoastPolygon::SetDownCoastThisIter(bool const bFlag)
{
   m_bDownCoastThisIter = bFlag;
}

bool CGeomCoastPolygon::bDownCoastThisIter(void) const
{
   return m_bDownCoastThisIter;
}


int CGeomCoastPolygon::nGetGlobalID(void) const
{
   return m_nGlobalID;
}

int CGeomCoastPolygon::nGetCoastID(void) const
{
   return m_nCoastID;
}

// void CGeomCoastPolygon::SetCoastNode(int const nNode)
// {
//    m_nCoastNode = nNode;
// }

int CGeomCoastPolygon::nGetNodeCoastPoint(void) const
{
   return m_nCoastNode;
}

CGeom2DIPoint* CGeomCoastPolygon::pPtiGetNode(void)
{
   return &m_PtiNode;

}

CGeom2DIPoint* CGeomCoastPolygon::pPtiGetAntiNode(void)
{
   return &m_PtiAntinode;
}

void CGeomCoastPolygon::SetNumCellsInPolygon(int const nCells)
{
   m_nNumCells = nCells;
}

int CGeomCoastPolygon::nGetNumCellsinPolygon(void) const
{
   return m_nNumCells;
}

int CGeomCoastPolygon::nGetUpCoastProfile(void) const
{
   return m_nProfileUpCoast;
}

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

CGeom2DPoint* CGeomCoastPolygon::pPtGetBoundaryPoint(int const nPoint)
{
   // NOTE no check to see if nPoint < m_VPoints.size()
   return &m_VPoints[nPoint];
}

int CGeomCoastPolygon::nGetBoundarySize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

int CGeomCoastPolygon::nGetUpCoastProfileNumPointsUsed(void) const
{
   return m_nProfileUpCoastNumPointsUsed;
}

int CGeomCoastPolygon::nGetDownCoastProfileNumPointsUsed(void) const
{
   return m_nProfileDownCoastNumPointsUsed;
}

void CGeomCoastPolygon::SetSeawaterVolume(const double dDepth)
{
   m_dSeawaterVolume = dDepth;
}

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

//! Re-initializes this timestep's eposition of unconsolidated fine sediment on this polygon
void CGeomCoastPolygon::SetZeroDepositionUnconsFine(void)
{
   m_dDepositionUnconsFine = 0;
}

//! Returns this timestep's deposition of fine unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetDepositionUnconsFine(void) const
{
   return m_dDepositionUnconsFine;
}

//! Re-initializes this timestep's deposition of unconsolidated sand sediment on this polygon
void CGeomCoastPolygon::SetZeroDepositionUnconsSand(void)
{
   m_dDepositionUnconsSand = 0;
}

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

//! Re-initializes this timestep's deposition of unconsolidated coarse sediment on this polygon
void CGeomCoastPolygon::SetZeroDepositionUnconsCoarse(void)
{
   m_dDepositionUnconsCoarse = 0;
}

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

void CGeomCoastPolygon::SetUpCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnUpCoastAdjacentPolygon = *pnVPolygons;
}

int CGeomCoastPolygon::nGetUpCoastAdjacentPolygon(int const nIndex) const
{
//    assert(nIndex < m_VnUpCoastAdjacentPolygon.size());
   return m_VnUpCoastAdjacentPolygon[nIndex];
}

int CGeomCoastPolygon::nGetNumUpCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnUpCoastAdjacentPolygon.size());
}


void CGeomCoastPolygon::SetDownCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnDownCoastAdjacentPolygon = *pnVPolygons;
}

int CGeomCoastPolygon::nGetDownCoastAdjacentPolygon(int const nIndex) const
{
//    assert(nIndex < m_VnDownCoastAdjacentPolygon.size());
   return m_VnDownCoastAdjacentPolygon[nIndex];
}

int CGeomCoastPolygon::nGetNumDownCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnDownCoastAdjacentPolygon.size());
}

void CGeomCoastPolygon::SetUpCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdUpCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

double CGeomCoastPolygon::dGetUpCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // NOTE no check to see if nIndex < m_VdUpCoastAdjacentPolygonBoundaryShare.size()
   return m_VdUpCoastAdjacentPolygonBoundaryShare[nIndex];
}

void CGeomCoastPolygon::SetDownCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdDownCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

double CGeomCoastPolygon::dGetDownCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // NOTE no check to see if nIndex < m_VdDownCoastAdjacentPolygonBoundaryShare.size()
   return m_VdDownCoastAdjacentPolygonBoundaryShare[nIndex];
}

int CGeomCoastPolygon::nGetPointInPolygonSearchStartPoint(void) const
{
   return m_nPointInPolygonSearchStartPoint;
}

void CGeomCoastPolygon::SetAvgUnconsD50(double const dD50)
{
   m_dAvgUnconsD50 = dD50;
}

double CGeomCoastPolygon::dGetAvgUnconsD50(void) const
{
   return m_dAvgUnconsD50;
}

void CGeomCoastPolygon::Display(void)
{
//    cout << endl;
//    for (int n = 0; n < static_cast<int>(m_VPoints.size()); n++)
//       cout << "[" << m_VPoints[n].dGetX() << "][" << m_VPoints[n].dGetY() << "], ";
//    cout << endl;
//    cout.flush();
}

void CGeomCoastPolygon::AddCircularity(int const nPoly)
{
   m_VnCircularityWith.push_back(nPoly);
}

vector<int> CGeomCoastPolygon::VnGetCircularities(void)
{
   return m_VnCircularityWith;
}

void CGeomCoastPolygon::AddCliffCollapseErosionFine(double const dDepth)
{
   m_dCliffCollapseErosionFine += dDepth;
}

double CGeomCoastPolygon::dGetCliffCollapseErosionFine(void) const
{
   return m_dCliffCollapseErosionFine;
}

void CGeomCoastPolygon::AddCliffCollapseErosionSand(double const dDepth)
{
   m_dCliffCollapseErosionSand += dDepth;
}

double CGeomCoastPolygon::dGetCliffCollapseErosionSand(void) const
{
   return m_dCliffCollapseErosionSand;
}

void CGeomCoastPolygon::AddCliffCollapseErosionCoarse(double const dDepth)
{
   m_dCliffCollapseErosionCoarse += dDepth;
}

double CGeomCoastPolygon::dGetCliffCollapseErosionCoarse(void) const   
{
   return m_dCliffCollapseErosionCoarse;
}

void CGeomCoastPolygon::AddCliffCollapseUnconsSandDeposition(double const dDepth)
{
   m_dCliffCollapseTalusSand += dDepth;
}

double CGeomCoastPolygon::dGetCliffCollapseUnconsSandDeposition(void) const
{
   return m_dCliffCollapseTalusSand;
}

void CGeomCoastPolygon::AddCliffCollapseUnconsCoarseDeposition(double const dDepth)
{  
   m_dCliffCollapseTalusCoarse += dDepth;
}

double CGeomCoastPolygon::dGetCliffCollapseUnconsCoarseDeposition(void) const
{
   return m_dCliffCollapseTalusCoarse;
}

void CGeomCoastPolygon::AddUnconsSandFromShorePlatform(double const dDepth)
{
   m_dSandFromPlatformErosion += dDepth;
}

double CGeomCoastPolygon::dGetUnconsSandFromShorePlatform(void) const
{
   return m_dSandFromPlatformErosion;
}

void CGeomCoastPolygon::AddUnconsCoarseFromShorePlatform(double const dDepth)
{
   m_dCoarseFromPlatformErosion += dDepth;
}

double CGeomCoastPolygon::dGetUnconsCoarseFromShorePlatform(void) const
{
   return m_dCoarseFromPlatformErosion;
}

void CGeomCoastPolygon::SetStoredUnconsFine(double const dDepth)
{
   m_dStoredUnconsFine = dDepth;
}

double CGeomCoastPolygon::dGetStoredUnconsFine(void) const
{
   return m_dStoredUnconsFine;
}

void CGeomCoastPolygon::SetStoredUnconsSand(double const dDepth)
{
   m_dStoredUnconsSand = dDepth;
}

double CGeomCoastPolygon::dGetStoredUnconsSand(void) const
{
   return m_dStoredUnconsSand;
}

void CGeomCoastPolygon::SetStoredUnconsCoarse(double const dDepth)
{
   m_dStoredUnconsCoarse = dDepth;
}

double CGeomCoastPolygon::dGetStoredUnconsCoarse(void) const
{
   return m_dStoredUnconsCoarse;
}

