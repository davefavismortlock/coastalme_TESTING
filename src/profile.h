/*!
 *
 * \class CGeomProfile
 * \brief Geometry class used to represent coast profile objects
 * \details TODO This is a more detailed description of the CGeomProfile class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file profile.h
 * \brief Contains CGeomProfile definitions
 *
 */

#ifndef PROFILE_H
#define PROFILE_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "2d_point.h"
#include "2di_point.h"
#include "multi_line.h"
#include "raster_grid.h"

class CGeomProfile : public CGeomMultiLine
{
private:
   bool
      m_bStartOfCoast,
      m_bEndOfCoast,
      m_bHitLand,
      m_bHitCoast,
      m_bTooShort,
      m_bTruncated,
      m_bHitAnotherProfile;
   int
      m_nNumCoastPoint;             // The coastline point at which this profile hits the coast (not necessarily coincident wih the profile start cell)
   double
      m_dDeepWaterWaveHeight,       // The wave height at the end of the profile
      m_dDeepWaterWaveAngle,  // The wave orientation at the end of the profile
      m_dDeepWaterWavePeriod;       // The wave period at the end of the profile

   vector<CGeom2DIPoint>
      m_VCellInProfile;             // In grid CRS, the integer coords of the cells 'under' this profile. NOTE Point zero is the same as 'cell marked as coastline' in coast object

   // The following have the same length as m_VCellInProfile
   vector<CGeom2DPoint>
      m_VCellInProfileExtCRS;       // In external CRS, the coords of cells 'under' this profile
//    vector<bool>
//       m_bVShared;                // Is this profile point part of a multi-line?
   vector<int>
      m_VnCoastPolyToLeft,          // The ID of the CoastPolygon to the left (looking seaward)
      m_VnCoastPolyToRight;         // The ID of the CoastPolygon to the right (looking seaward)

public:
   explicit CGeomProfile(int const);
   ~CGeomProfile(void);

   int nGetNumCoastPoint(void) const;

   void SetStartOfCoast(bool const);
   bool bStartOfCoast(void) const;

   void SetEndOfCoast(bool const);
   bool bEndOfCoast(void) const;

   void SetHitLand(bool const);
   bool bHitLand(void) const;

   void SetHitCoast(bool const);
   bool bHitCoast(void) const;

   void SetTooShort(bool const);
   bool bTooShort(void) const;

   void SetTruncated(bool const);
   bool bTruncated(void) const;

   void SetHitAnotherProfile(bool const);
   bool bHitAnotherProfile(void) const;

   bool bProfileOK(void) const;
   bool bOKIncStartAndEndOfCoast(void) const;
   bool bOKIncStartOfCoast(void) const;

   void SetAllPointsInProfile(vector<CGeom2DPoint> const*);
   void SetPointInProfile(int const, double const, double const);
   void AppendPointInProfile(double const, double const);
   void AppendPointInProfile(CGeom2DPoint const*);
   void TruncateProfile(int const);
//    void TruncateAndSetPointInProfile(int const, double const, double const);
   bool bInsertIntersection(double const, double const, int const);
//    void ShowProfile(void) const;
   int nGetProfileSize(void) const;
   CGeom2DPoint* pPtGetPointInProfile(int const);
   vector<CGeom2DPoint> PtVGetThisPointAndAllAfter(int const);
   void RemoveLineSegment(int const);
   bool bIsPointInProfile(double const, double const);
   bool bIsPointInProfile(double const, double const, int&);
//    int nFindInsertionLineSeg(double const, double const);

//    void AppendPointShared(bool const);
//    bool bPointShared(int const) const;

//    void SetCoastPolyToLeft(int const, int const);
//    int nGetCoastPolyToleft(int const);

//    void SetCoastPolyToRight(int const, int const);
//    int nGetCoastPolyToRight(int const);

   void AppendCellInProfile(CGeom2DIPoint const*);
   void AppendCellInProfile(int const, int const);
//    void SetCellsInProfile(vector<CGeom2DIPoint>*);
   vector<CGeom2DIPoint>* pPtiVGetCellsInProfile(void);
   CGeom2DIPoint* pPtiGetCellInProfile(int const);
   int nGetNumCellsInProfile(void) const;

   void AppendCellInProfileExtCRS(double const, double const);
   void AppendCellInProfileExtCRS(CGeom2DPoint const*);
//    vector<CGeom2DPoint>* PtVGetCellsInProfileExtCRS(void);

   int nGetCellGivenDepth(CGeomRasterGrid const*, double const);

   void SetProfileDeepWaterWaveHeight(double const);
   double dGetProfileDeepWaterWaveHeight(void) const;

   void SetProfileDeepWaterWaveAngle(double const);
   double dGetProfileDeepWaterWaveAngle(void) const;

   void SetProfileDeepWaterWavePeriod(double const);
   double dGetProfileDeepWaterWavePeriod(void) const;
};
#endif //PROFILE_H

