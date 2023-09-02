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

//! Constructor
CRWCliff::CRWCliff(CRWCoast* pCoastIn, int const nCoast, int const nPointOnCoast, double const dCellSide, double const dNotchDepthIn, double const dNotchElevIn, double const dAccumWaveEnergyIn)
{
   m_bCliffCollapse      =
   m_bAllSedimentGone    = false;

   pCoast                = pCoastIn;

   m_nCoast              = nCoast;
   m_nPointOnCoast       = nPointOnCoast;
   m_nCategory           = LF_CAT_CLIFF;

   m_dMaxDepth           = dCellSide;
   m_dNotchDepth         = dNotchDepthIn;
   m_dNotchBaseElev      = dNotchElevIn;
   m_dTotAccumWaveEnergy = dAccumWaveEnergyIn;

//    assert(m_dRemaining >=0);
}

//! Destructor
CRWCliff::~CRWCliff(void)
{
}

// bool CRWCliff::bHasCollapsed(void) const
// {
//    return m_bCliffCollapse;
// }

//! Sets the cliff collapse switch
void CRWCliff::SetCliffCollapse(bool const bStatus)
{
   m_bCliffCollapse = bStatus;
}

//! Returns the value of the all sediment gone switch
bool CRWCliff::bAllSedimentGone(void) const
{
   return m_bAllSedimentGone;
}

//! Sets the all sediment gone switch to true
void CRWCliff::SetAllSedimentGone(void)
{
   m_bAllSedimentGone = true;
}

//! Returns the elevation of the base of the erosional notch
double CRWCliff::dGetNotchBaseElev(void) const
{
   return m_dNotchBaseElev;
}

//! Sets the elevation of the base of the erosional notch
void CRWCliff::SetNotchBaseElev(double const dNewElev)
{
   m_dNotchBaseElev = dNewElev;
}

// void CRWCliff::SetRemaining(double const dLenIn)
// {
//    m_dRemaining = dLenIn;
// }

//! Returns the XY-plane length (in external CRS units) of the remaining sediment on the coast cell at the elevation of the notch. Note that notch depth is not changed
double CRWCliff::dGetRemaining(void) const
{
   return (m_dMaxDepth - m_dNotchDepth);
}

//! Sets the depth of the erosional notch. Note that the remaining sediment on the coast cell at the elevation of the notch is not changed
void CRWCliff::SetNotchDepth(double const dLenIn)
{
   m_dNotchDepth = dLenIn;
}

//! Returns the depth of the erosional notch
double CRWCliff::dGetNotchDepth(void) const
{
   return m_dNotchDepth;
}

//! Returns true if the notch has reached the edge of the cell, or if the notch overhang exceeds the critical notch overhang
bool CRWCliff::bReadyToCollapse(double const dThresholdNotchDepth) const
{
   if (m_dNotchDepth >= dThresholdNotchDepth)
      return true;
   else
      return false;
}

//! Increases the XY-plane length (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
void CRWCliff::DeepenErosionalNotch(double const dLenIn)
{
   m_dNotchDepth += dLenIn;

   // Constrain the notch depth, it cannot be greater than the max notch depth
   m_dNotchDepth = tMin(m_dNotchDepth, m_dMaxDepth);
   
//    assert((m_dMaxDepth - m_dNotchDepth) >=0);
}

//! Display (dummy)
void CRWCliff::Display(void)
{
   cout << endl;
//    for (int n = 0; n < static_cast<int>(m_VPoints.size()); n++)
//       cout << "[" << m_VPoints[n].dGetX() << "][" << m_VPoints[n].dGetY() << "], ";
//    cout << endl;
//    cout.flush();
}
