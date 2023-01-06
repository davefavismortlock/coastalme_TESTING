/*!
 *
 * \file coast.cpp
 * \brief CRWCoast routines
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2021
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

#include <vector>
#include <algorithm>

#include "cme.h"
#include "coast.h"
#include "line.h"
#include "i_line.h"

CRWCoast::CRWCoast(void)
    : m_nSeaHandedness(NULL_HANDED),
      m_nStartEdge(INT_NODATA),
      m_nEndEdge(INT_NODATA),
      m_dCurvatureDetailedMean(0),
      m_dCurvatureDetailedSTD(0),
      m_dCurvatureSmoothMean(0),
      m_dCurvatureSmoothSTD(0)
{
}

CRWCoast::~CRWCoast(void)
{
   for (unsigned int i = 0; i < m_pVLandforms.size(); i++)
      delete m_pVLandforms[i];

   for (unsigned int i = 0; i < m_pVPolygon.size(); i++)
      delete m_pVPolygon[i];
}

void CRWCoast::SetSeaHandedness(int const nNewHandedness)
{
   m_nSeaHandedness = nNewHandedness;
}

int CRWCoast::nGetSeaHandedness(void) const
{
   return m_nSeaHandedness;
}

void CRWCoast::SetStartEdge(int const nEdge)
{
   m_nStartEdge = nEdge;
}

int CRWCoast::nGetStartEdge(void) const
{
   return m_nStartEdge;
}

void CRWCoast::SetEndEdge(int const nEdge)
{
   m_nEndEdge = nEdge;
}

int CRWCoast::nGetEndEdge(void) const
{
   return m_nEndEdge;
}

void CRWCoast::SetCoastlineExtCRS(CGeomLine const *pLCoast)
{
   m_LCoastlineExtCRS = *pLCoast;

   int nLen = m_LCoastlineExtCRS.nGetSize();

   m_VnProfileNumber = vector<int>(nLen, INT_NODATA);
   m_VnPolygonNode = vector<int>(nLen, INT_NODATA);
   m_VnBreakingDistance = vector<int>(nLen, INT_NODATA);

   m_VdCurvatureDetailed = vector<double>(nLen, DBL_NODATA);
   m_VdCurvatureSmooth = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWaveAngle = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWavePeriod = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdWaveSetupSurge = vector<double>(nLen, 0); // it is better to initiate with DBL_NODATA but some values are outside of range in the interpolation
   // m_VdStormSurge = vector<double>(nLen, DBL_NODATA);
   m_VdRunUp = vector<double>(nLen, 0);
   m_VdCoastWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveAngle = vector<double>(nLen, DBL_NODATA);
   m_VdDepthOfBreaking = vector<double>(nLen, DBL_NODATA);
   m_VdFluxOrientation = vector<double>(nLen, DBL_NODATA);
   m_VdWaveEnergyAtBreaking = vector<double>(nLen, 0);
}

void CRWCoast::AppendPointToCoastlineExtCRS(double const dX, double const dY)
{
   // Appends a coastline point (in external CRS), also appends dummy values for curvature, breaking wave height, wave angle, and flux orientation
   m_LCoastlineExtCRS.Append(dX, dY);

   m_VnProfileNumber.push_back(INT_NODATA);
   m_VnPolygonNode.push_back(INT_NODATA);
   m_VnBreakingDistance.push_back(INT_NODATA);

   m_VdCurvatureDetailed.push_back(DBL_NODATA);
   m_VdCurvatureSmooth.push_back(DBL_NODATA);
   m_VdDeepWaterWaveHeight.push_back(DBL_NODATA);
   m_VdDeepWaterWaveAngle.push_back(DBL_NODATA);
   m_VdBreakingWaveHeight.push_back(DBL_NODATA);
   m_VdWaveSetupSurge.push_back(DBL_NODATA);
   // m_VdStormSurge.push_back(DBL_NODATA);
   m_VdRunUp.push_back(DBL_NODATA);
   m_VdBreakingWaveAngle.push_back(DBL_NODATA);
   m_VdDepthOfBreaking.push_back(DBL_NODATA);
   m_VdFluxOrientation.push_back(DBL_NODATA);
   m_VdWaveEnergyAtBreaking.push_back(0);
}

// void CRWCoast::SetFloodWaveSetupPointExtCRS(CGeomLine const *pLCoast)
// {
//    m_LFloodWaveSetupExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupLineExtCRS.nGetSize();
// }

// void CRWCoast::SetFloodWaveSetupSurgePointExtCRS(CGeomLine const *pLCoast)
// {
//    m_LFloodWaveSetupSurgeExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupSurgeLineExtCRS.nGetSize();
// }

// void CRWCoast::SetFloodWaveSetupSurgeRunupPointExtCRS(CGeomLine const *pLCoast)
// {
//    m_LFloodWaveSetupSurgeRunupExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupSurgeRunupLineExtCRS.nGetSize();
// }

CGeomLine *CRWCoast::pLGetCoastlineExtCRS(void)
{
   return &m_LCoastlineExtCRS;
}

CGeom2DPoint *CRWCoast::pPtGetCoastlinePointExtCRS(int const n)
{
   // Point is in external CRS NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
   return &m_LCoastlineExtCRS[n];
}

// CGeomLine *CRWCoast::pLGetFloodWaveSetupExtCRS(void)
// {
//    return &m_LFloodWaveSetupExtCRS;
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupPointExtCRS(int const n)
// {
//    // Point is in external CRS NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupExtCRS[n];
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupSurgePointExtCRS(int const n)
// {
//    // Point is in external CRS NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupSurgeExtCRS[n];
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupSurgeRunupPointExtCRS(int const n)
// {
//    // Point is in external CRS NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupSurgeRunupExtCRS[n];
// }

int CRWCoast::nGetCoastlineSize(void) const
{
   return m_LCoastlineExtCRS.nGetSize();
}

// void CRWCoast::DisplayCoastline(void)
// {
//    m_LCoastlineExtCRS.Display();
// }

void CRWCoast::SetCoastlineGridCRS(CGeomILine const *pILCoastCells)
{
   m_ILCellsMarkedAsCoastline = *pILCoastCells;
}

// void CRWCoast::AppendCellMarkedAsCoastline(CGeom2DIPoint const* pPti)
// {
//    m_ILCellsMarkedAsCoastline.Append(*pPti);
// }
//
// void CRWCoast::AppendCellMarkedAsCoastline(int const nX, int const nY)
// {
//    m_ILCellsMarkedAsCoastline.Append(CGeom2DIPoint(nX, nY));
// }

CGeom2DIPoint *CRWCoast::pPtiGetCellMarkedAsCoastline(int const n)
{
   // NOTE No check to see if n < size()
   return &m_ILCellsMarkedAsCoastline[n];
}

// int CRWCoast::nGetNCellsMarkedAsCoastline(void) const
// {
//    return m_ILCellsMarkedAsCoastline.size();
// }

// double CRWCoast::dGetCoastlineSegmentLength(int const m, int const n)
// {
//    // NOTE no check to see that m is < m_LCoastlineExtCRS.Size(), same for n
//    if (m == n)
//       return 0;
//
//    return hypot(m_LCoastlineExtCRS[n].dGetX() - m_LCoastlineExtCRS[m].dGetX(), m_LCoastlineExtCRS[n].dGetY() - m_LCoastlineExtCRS[m].dGetY());
// }

// double CRWCoast::dGetCoastlineLengthSoFar(int const n)
// {
//    // NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
//    double dLen = 0;
//    for (int m = 0; m < n; m++)
//       dLen += dGetCoastlineSegmentLength(m, m+1);
//    return dLen;
// }

//! Returns the coastline number given a cell, or INT_NODATA if neither this cell or any of its neighbouring cells are 'under' a coastline. If it is a neighbouring cell that is under the coastline, then it also changes the cell that is supplied as an input parameter
int CRWCoast::nGetCoastPointGivenCell(CGeom2DIPoint *pPtiCell)
{
   for (int nCoastPoint = 0; nCoastPoint < m_ILCellsMarkedAsCoastline.nGetSize(); nCoastPoint++)
   {
      if (m_ILCellsMarkedAsCoastline[nCoastPoint] == pPtiCell)
      {
         return nCoastPoint;
      }
   }

   // This cell is not under a coastline, so try the adjacent cells
   int
       n = -1,
       nX = pPtiCell->nGetX(),
       nY = pPtiCell->nGetY(),
       nXAdj = 0,
       nYAdj = 0;

   while (n <= 7)
   {
      switch (++n)
      {
      case 0:
         nXAdj = nX;
         nYAdj = nY - 1;
         break;
      case 1:
         nXAdj = nX + 1;
         nYAdj = nY - 1;
         break;
      case 2:
         nXAdj = nX + 1;
         nYAdj = nY;
         break;
      case 3:
         nXAdj = nX + 1;
         nYAdj = nY + 1;
         break;
      case 4:
         nXAdj = nX;
         nYAdj = nY + 1;
         break;
      case 5:
         nXAdj = nX - 1;
         nYAdj = nY + 1;
         break;
      case 6:
         nXAdj = nX - 1;
         nYAdj = nY;
         break;
      case 7:
         nXAdj = nX - 1;
         nYAdj = nY - 1;
         break;
      }

      CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
      for (int nCoastPoint = 0; nCoastPoint < m_ILCellsMarkedAsCoastline.nGetSize(); nCoastPoint++)
      {
         if (m_ILCellsMarkedAsCoastline[nCoastPoint] == &PtiTmp)
         {
            *pPtiCell = PtiTmp;
            return nCoastPoint;
         }
      }
   }

   return INT_NODATA;
}

double CRWCoast::dGetDetailedCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VdCurvatureDetailed.Size()
   return m_VdCurvatureDetailed[nCoastPoint];
}

void CRWCoast::SetDetailedCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_VdCurvatureDetailed.size()
   m_VdCurvatureDetailed[nCoastPoint] = dCurvature;
}

vector<double> *CRWCoast::pVGetDetailedCurvature(void)
{
   return &m_VdCurvatureDetailed;
}

double CRWCoast::dGetSmoothCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VdCurvatureSmooth.Size()
   return m_VdCurvatureSmooth[nCoastPoint];
}

void CRWCoast::SetSmoothCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_VdCurvatureSmooth.size()
   m_VdCurvatureSmooth[nCoastPoint] = dCurvature;
}

vector<double> *CRWCoast::pVGetSmoothCurvature(void)
{
   return &m_VdCurvatureSmooth;
}

void CRWCoast::SetDetailedCurvatureMean(double const dMean)
{
   m_dCurvatureDetailedMean = dMean;
}

double CRWCoast::dGetDetailedCurvatureMean(void) const
{
   return m_dCurvatureDetailedMean;
}

void CRWCoast::SetDetailedCurvatureSTD(double const dSTD)
{
   m_dCurvatureDetailedSTD = dSTD;
}

double CRWCoast::dGetDetailedCurvatureSTD(void) const
{
   return m_dCurvatureDetailedSTD;
}

void CRWCoast::SetSmoothCurvatureMean(double const dMean)
{
   m_dCurvatureSmoothMean = dMean;
}

double CRWCoast::dGetSmoothCurvatureMean(void) const
{
   return m_dCurvatureSmoothMean;
}

void CRWCoast::SetSmoothCurvatureSTD(double const dSTD)
{
   m_dCurvatureSmoothSTD = dSTD;
}

double CRWCoast::dGetSmoothCurvatureSTD(void) const
{
   return m_dCurvatureSmoothSTD;
}

CGeomProfile *CRWCoast::pGetProfile(int const nProfile)
{
   // NOTE No safety check that nProfile < m_VProfile.size()
   return &m_VProfile[nProfile];
}

void CRWCoast::AppendProfile(int const nCoastPoint, int const nProfile)
{
   CGeomProfile Profile(nCoastPoint);
   m_VProfile.push_back(Profile);

   m_VnProfileNumber[nCoastPoint] = nProfile;
}

// void CRWCoast::ReplaceProfile(int const nProfile, vector<CGeom2DPoint> const* pPtVProfileNew)
// {
//    // NOTE No safety check that nProfile < m_VProfile.size()
//    m_VProfile[nProfile].SetAllPointsInProfile(pPtVProfileNew);
// }

int CRWCoast::nGetNumProfiles(void) const
{
   return static_cast<int>(m_VProfile.size());
}

bool CRWCoast::bIsProfileStartPoint(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VnProfileNumber.Size()
   if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
      return true;

   return false;
}

int CRWCoast::nGetProfileNumber(int const nCoastPoint) const
{
   // Returns INT_NODATA if no profile at this point
   return m_VnProfileNumber[nCoastPoint];
}

void CRWCoast::CreateAlongCoastProfileIndex(void)
{
   // Creates an index containing the numbers of the coastline-normal profiles in along-coast sequence
   for (int nCoastPoint = 0; nCoastPoint < m_LCoastlineExtCRS.nGetSize(); nCoastPoint++)
   {
      if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
         m_VnProfileCoastIndex.push_back(m_VnProfileNumber[nCoastPoint]);
   }
}

int CRWCoast::nGetProfileFromAlongCoastProfileIndex(int const n) const
{
   // Returns the number of the coastline-normal profile which is at position n in the along-coast sequence
   return m_VnProfileCoastIndex[n];
}

int CRWCoast::nGetDownCoastProfileNumber(int const nProfile) const
{
   // Return the number of the profile which is adjacent to and down-coast from the specified profile. It returns INT_NODATA if there is no valid up-coast profile
   for (unsigned int n = 0; n < m_VnProfileCoastIndex.size() - 1; n++)
   {
      if (nProfile == m_VnProfileCoastIndex[n])
         return m_VnProfileCoastIndex[n + 1];
   }

   // At end of m_VnProfileCoastIndex, so no down-coast profile
   return INT_NODATA;
}

// int CRWCoast::nGetAlongCoastlineIndexOfProfile(int const nProfile)
// {
//    // Returns the along-coastline index of a coastline-normal profile
//    for (int n = 0; n < m_VnProfileCoastIndex.size(); n++)
//       if (m_VnProfileCoastIndex[n] == nProfile)
//          return n;
//    return -1;
// }

void CRWCoast::SetCoastDeepWaterWaveHeight(int const nCoastPoint, double const dHeight)
{
   // NOTE no check to see if nCoastPoint < m_VdDeepWaterWaveHeight.size()
   m_VdDeepWaterWaveHeight[nCoastPoint] = dHeight;
}

double CRWCoast::dGetCoastDeepWaterWaveHeight(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdDeepWaterWaveHeight.size()
   return m_VdDeepWaterWaveHeight[nCoastPoint];
}

void CRWCoast::SetCoastDeepWaterWaveAngle(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_VdDeepWaterWaveAngle.size()
   m_VdDeepWaterWaveAngle[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetCoastDeepWaterWaveAngle(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdDeepWaterWaveAngle.size()
   return m_VdDeepWaterWaveAngle[nCoastPoint];
}

void CRWCoast::SetCoastDeepWaterWavePeriod(int const nCoastPoint, double const dPeriod)
{
   m_VdDeepWaterWavePeriod[nCoastPoint] = dPeriod;
}

double CRWCoast::dGetCoastDeepWaterWavePeriod(int const nCoastPoint) const
{
   return m_VdDeepWaterWavePeriod[nCoastPoint];
}

void CRWCoast::SetBreakingWaveHeight(int const nCoastPoint, double const dHeight)
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   m_VdBreakingWaveHeight[nCoastPoint] = dHeight;
}

double CRWCoast::dGetBreakingWaveHeight(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   return m_VdBreakingWaveHeight[nCoastPoint];
}

void CRWCoast::SetWaveSetupSurge(int const nCoastPoint, double const dWaveSetup)
{
   m_VdWaveSetupSurge[nCoastPoint] = dWaveSetup;
}

double CRWCoast::dGetWaveSetupSurge(int const nCoastPoint) const
{
   return m_VdWaveSetupSurge[nCoastPoint];
}

// void CRWCoast::SetStormSurge(int const nCoastPoint, double const dStormSurge)
// {
//    m_VdStormSurge[nCoastPoint] = dStormSurge;
// }

// double CRWCoast::dGetStormSurge(int const nCoastPoint) const
// {
//    return m_VdStormSurge[nCoastPoint];
// }

void CRWCoast::SetRunUp(int const nCoastPoint, double const dRunUp)
{
   m_VdRunUp[nCoastPoint] = dRunUp;
}

double CRWCoast::dGetRunUp(int const nCoastPoint) const
{
   return m_VdRunUp[nCoastPoint];
}

double CRWCoast::dGetLevel(int const nCoastPoint, int const level) const
{
   switch (level)
   {
   case 0: // WAVESETUPSURGE:
      return m_VdWaveSetupSurge[nCoastPoint];
      break;
   case 1: // WAVESETUPSURGE + RUNUP:
      return m_VdWaveSetupSurge[nCoastPoint] + m_VdRunUp[nCoastPoint];
      break;
   default:
      return 0;
   }
}

void CRWCoast::SetCoastWaveHeight(int const nCoastPoint, double const dHeight)
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   m_VdCoastWaveHeight[nCoastPoint] = dHeight;
}

double CRWCoast::dGetCoastWaveHeight(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   return m_VdCoastWaveHeight[nCoastPoint];
}

void CRWCoast::SetBreakingWaveAngle(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   m_VdBreakingWaveAngle[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetBreakingWaveAngle(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   return m_VdBreakingWaveAngle[nCoastPoint];
}

void CRWCoast::SetDepthOfBreaking(int const nCoastPoint, double const dDepth)
{
   // NOTE no check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   m_VdDepthOfBreaking[nCoastPoint] = dDepth;
}

double CRWCoast::dGetDepthOfBreaking(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   return m_VdDepthOfBreaking[nCoastPoint];
}

void CRWCoast::SetBreakingDistance(int const nCoastPoint, int const nDist)
{
   // NOTE no check to see if nCoastPoint < m_VnBreakingDistance.size()
   m_VnBreakingDistance[nCoastPoint] = nDist;
}

int CRWCoast::nGetBreakingDistance(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VnBreakingDistance.size()
   return m_VnBreakingDistance[nCoastPoint];
}

void CRWCoast::SetFluxOrientation(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_VdFluxOrientation.size()
   m_VdFluxOrientation[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetFluxOrientation(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdFluxOrientation.size()
   return m_VdFluxOrientation[nCoastPoint];
}

void CRWCoast::SetWaveEnergyAtBreaking(int const nCoastPoint, double const dEnergy)
{
   // NOTE no check to see if nCoastPoint < m_VdWaveEnergyAtBreaking.size()
   //    assert(isfinite(dEnergy));
   m_VdWaveEnergyAtBreaking[nCoastPoint] = dEnergy;
}

double CRWCoast::dGetWaveEnergyatBreaking(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdWaveEnergyAtBreaking.size()
   //    assert(isfinite(m_VdWaveEnergyAtBreaking[nCoastPoint]));
   return m_VdWaveEnergyAtBreaking[nCoastPoint];
}

void CRWCoast::AppendCoastLandform(CACoastLandform *pCoastLandform)
{
   m_pVLandforms.push_back(pCoastLandform);
}

CACoastLandform *CRWCoast::pGetCoastLandform(int const nCoastPoint)
{
   if (nCoastPoint < static_cast<int>(m_pVLandforms.size()))
      return m_pVLandforms[nCoastPoint];

   return NULL;
}

void CRWCoast::SetPolygonNode(int const nPoint, int const nNode)
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   m_VnPolygonNode[nPoint] = nNode;
}

int CRWCoast::nGetPolygonNode(int const nPoint) const
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   return m_VnPolygonNode[nPoint];
}

void CRWCoast::CreatePolygon(int const nGlobalID, int const nCoastID, int const nCoastPoint, CGeom2DIPoint const *PtiNode, CGeom2DIPoint const *PtiAntiNode, int const nProfileUpCoast, int const nProfileDownCoast, vector<CGeom2DPoint> const *pVIn, int const nPointsUpCoastProfile, int const nPointsDownCoastProfile, int const nPointInPolygonStartPoint)
{
   CGeomCoastPolygon *pPolygon = new CGeomCoastPolygon(nGlobalID, nCoastID, nCoastPoint, nProfileUpCoast, nProfileDownCoast, pVIn, nPointsUpCoastProfile, nPointsDownCoastProfile, PtiNode, PtiAntiNode, nPointInPolygonStartPoint);

   m_pVPolygon.push_back(pPolygon);
}

int CRWCoast::nGetNumPolygons(void) const
{
   return static_cast<int>(m_pVPolygon.size());
}

CGeomCoastPolygon *CRWCoast::pGetPolygon(int const nPoly) const
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   return m_pVPolygon[nPoly];
}

void CRWCoast::AppendPolygonLength(const double dLength)
{
   m_VdPolygonLength.push_back(dLength);
}

double CRWCoast::dGetPolygonLength(int const nIndex) const
{
   // NOTE no check to see if nIndex < m_VdPolygonLength.size()
   return m_VdPolygonLength[nIndex];
}

int CRWCoast::nGetNumShadowBoundaries(void)
{
   return static_cast<int>(m_LShadowBoundary.size());
}

void CRWCoast::AppendShadowBoundary(CGeomLine const *pLBoundary)
{
   m_LShadowBoundary.push_back(*pLBoundary);
}

CGeomLine *CRWCoast::pGetShadowBoundary(int const n)
{
   // NOTE no check to see if n < m_LShadowBoundary.size()
   return &m_LShadowBoundary[n];
}

int CRWCoast::nGetNumShadowDowndriftBoundaries(void)
{
   return static_cast<int>(m_LShadowDowndriftBoundary.size());
}

void CRWCoast::AppendShadowDowndriftBoundary(CGeomLine const *pLBoundary)
{
   m_LShadowDowndriftBoundary.push_back(*pLBoundary);
}

CGeomLine *CRWCoast::pGetShadowDowndriftBoundary(int const n)
{
   // NOTE no check to see if n < m_LShadowDowndriftBoundary.size()
   return &m_LShadowDowndriftBoundary[n];
}
