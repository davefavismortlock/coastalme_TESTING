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
// #include <assert.h>

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
//    m_nNumCells(0),
   m_nPointInPolygonSearchStartPoint(nPointInPolygonStartPoint),
//    m_dSeawaterVolume(0),
   m_dAvgUnconsD50(0),
   m_dDeltaPotentialTotalSediment(0),
   m_dDeltaEstimatedUnconsFine(0),
   m_dDeltaEstimatedUnconsSand(0),
   m_dDeltaEstimatedUnconsCoarse(0),
   m_dDeltaActualUnconsFine(0),
   m_dDeltaActualUnconsSand(0),
   m_dDeltaActualUnconsCoarse(0),
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


// void CGeomCoastPolygon::SetNumCellsInPolygon(int const nCells)
// {
//    m_nNumCells = nCells;
// }
//
//
// int CGeomCoastPolygon::nGetNumCellsinPolygon(void) const
// {
//    return m_nNumCells;
// }


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


// void CGeomCoastPolygon::SetSeawaterVolume(const double dDepth)
// {
//    m_dSeawaterVolume = dDepth;
// }

// double CGeomCoastPolygon::dGetSeawaterVolume(void) const
// {
//    return m_dSeawaterVolume;
// }

//! Adds a change in potential erosion or deposition to this timestep's total change in depth of unconsolidated sediment (all size classes) on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::AddDeltaPotentialTotalSediment(double const dDepth)
{
   m_dDeltaPotentialTotalSediment += dDepth;
}

//! Returns this timestep's total change in depth of unconsolidated sediment (all size classes) on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaPotentialErosion(void) const
{
   return m_dDeltaPotentialTotalSediment;
}

//! Sets a value for this timestep's estimated total change in depth of fine unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::SetDeltaEstimatedUnconsFine(double const dDepth)
{
   m_dDeltaEstimatedUnconsFine = dDepth;
}

//! Returns this timestep's estimate of total change in depth of fine unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaEstimatedUnconsFine(void) const
{
   return m_dDeltaEstimatedUnconsFine;
}

//! Sets a value for this timestep's estimated total change in depth of sand-sized unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::SetDeltaEstimatedUnconsSand(double const dDepth)
{
   m_dDeltaEstimatedUnconsSand = dDepth;
}

//! Returns this timestep's estimate of total change in depth of sand-sized unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaEstimatedUnconsSand(void) const
{
   return m_dDeltaEstimatedUnconsSand;
}

//! Sets a value for this timestep's estimate of total change in depth of coarse unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::SetDeltaEstimatedUnconsCoarse(double const dDepth)
{
   m_dDeltaEstimatedUnconsCoarse = dDepth;
}

//! Returns this timestep's estimate of total change in depth of coarse unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaEstimatedUnconsCoarse(void) const
{
   return m_dDeltaEstimatedUnconsCoarse;
}

//! Adds a change in erosion or deposition to this timestep's total actual change in depth of unconsolidated fine sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::AddDeltaActualUnconsFine(double const dDepth)
{
   m_dDeltaActualUnconsFine += dDepth;
}

// void CGeomCoastPolygon::SetDeltaActualUnconsFine(double const dDepth)
// {
//    m_dDeltaActualUnconsFine = dDepth;
// }

//! Returns this timestep's actual total change in depth of fine unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaActualUnconsFine(void) const
{
   return m_dDeltaActualUnconsFine;
}

//! Adds a change in erosion or deposition to this timestep's total actual change in depth of unconsolidated sand-sized sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::AddDeltaActualUnconsSand(double const dDepth)
{
   m_dDeltaActualUnconsSand += dDepth;
}

// void CGeomCoastPolygon::SetDeltaActualUnconsSand(double const dDepth)
// {
//    m_dDeltaActualUnconsSand = dDepth;
// }

//! Returns this timestep's actual total change in depth of sand-sized unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaActualUnconsSand(void) const
{
   return m_dDeltaActualUnconsSand;
}

//! Adds a change in erosion or deposition to this timestep's total actual change in depth of unconsolidated coarse sediment on this polygon (-ve erosion, +ve deposition)
void CGeomCoastPolygon::AddDeltaActualUnconsCoarse(double const dDepth)
{
   m_dDeltaActualUnconsCoarse += dDepth;
}

// void CGeomCoastPolygon::SetDeltaActualUnconsCoarse(double const dDepth)
// {
//    m_dDeltaActualUnconsCoarse = dDepth;
// }

//! Returns this timestep's actual total change in depth of coarse unconsolidated sediment on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaActualUnconsCoarse(void) const
{
   return m_dDeltaActualUnconsCoarse;
}

//! Returns this timestep's actual total change in depth of unconsolidated sediment (all size classes) on this polygon (-ve erosion, +ve deposition)
double CGeomCoastPolygon::dGetDeltaActualTotalSediment(void) const
{
   return m_dDeltaActualUnconsFine + m_dDeltaActualUnconsSand + m_dDeltaActualUnconsCoarse;
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
