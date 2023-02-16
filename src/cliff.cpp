/*!
 *
 * \file cliff.cpp
 * \brief CRWCliff routines
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
//#include <assert.h>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include "cme.h"
#include "cliff.h"


CRWCliff::CRWCliff(CRWCoast* pCoastIn, int const nCoast, int const nPointOnCoast, double const dRemainIn, double const dNotchElevIn, double const dNotchOverhangIn, double const dAccumWaveEnergyIn)
{
   m_bCliffCollapse   =
   m_bAllSedimentGone = false;

   pCoast = pCoastIn;

   m_nCoast        = nCoast;
   m_nPointOnCoast = nPointOnCoast;
   m_nCategory     = LF_CAT_CLIFF;

   m_dRemaining          = dRemainIn;
   m_dNotchBaseElev      = dNotchElevIn;
   m_dNotchOverhang      = dNotchOverhangIn;
   m_dTotAccumWaveEnergy = dAccumWaveEnergyIn;

//    assert(m_dRemaining >=0);
}

CRWCliff::~CRWCliff(void)
{
}

// bool CRWCliff::bHasCollapsed(void) const
// {
//    return m_bCliffCollapse;
// }

void CRWCliff::SetCliffCollapse(bool const bStatus)
{
   m_bCliffCollapse = bStatus;
}

bool CRWCliff::bAllSedimentGone(void) const
{
   return m_bAllSedimentGone;
}

void CRWCliff::SetAllSedimentGone(void)
{
   m_bAllSedimentGone = true;
}

double CRWCliff::dGetNotchBaseElev(void) const
{
   return m_dNotchBaseElev;
}

void CRWCliff::SetNotchBaseElev(double const dNewElev)
{
   m_dNotchBaseElev = dNewElev;
}

// void CRWCliff::SetRemaining(double const dLenIn)
// {
//    m_dRemaining = dLenIn;
// }

double CRWCliff::dGetRemaining(void) const
{
   return m_dRemaining;
}

void CRWCliff::SetNotchOverhang(double const dLenIn)
{
   m_dNotchOverhang = dLenIn;
}

double CRWCliff::dGetNotchOverhang(void) const
{
   return m_dNotchOverhang;
}

//! Returns true if the notch has reached the edge of the cell, or if the notch overhang exceeds the critical notch overhang
bool CRWCliff::bReadyToCollapse(double const dThresholdOverhang) const
{
   if ((m_dRemaining <= 0) || (m_dNotchOverhang >= dThresholdOverhang))
      return true;
   else
      return false;
}

// Increase the XY-plane length (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
double CRWCliff::dErodeNotch(double const dLenIn)
{
   // First constrain the supplied XY-plane length increment, it cannot exceed the XY-plane length of sediment remaining
   double dToRemove = tMin(m_dRemaining, dLenIn);
   m_dRemaining -= dToRemove;
   m_dNotchOverhang += dToRemove;

//    assert(m_dRemaining >=0);

   // Return the (possibly reduced) XY-plane length increment
   return dToRemove;
}

void CRWCliff::Display(void)
{
   cout << endl;
//    for (int n = 0; n < static_cast<int>(m_VPoints.size()); n++)
//       cout << "[" << m_VPoints[n].dGetX() << "][" << m_VPoints[n].dGetY() << "], ";
//    cout << endl;
//    cout.flush();
}
