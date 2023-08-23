/*!
 *
 * \file cell_layer.cpp
 * \brief CRWCellLayer routines
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
#include "cell_layer.h"

//! Constructor
CRWCellLayer::CRWCellLayer(void)
{
}

//! Destructor
CRWCellSediment* CRWCellLayer::pGetUnconsolidatedSediment(void)
{
   return &m_UnconsolidatedSediment;
}

CRWCellSediment* CRWCellLayer::pGetConsolidatedSediment(void)
{
   return &m_ConsolidatedSediment;
}

void CRWCellLayer::RemoveCliff(void)
{
   double
      dNow = 0,
      dLost = 0;

   dNow  = m_UnconsolidatedSediment.dGetFine();
   dLost = m_UnconsolidatedSediment.dGetNotchFineLost();
   m_UnconsolidatedSediment.SetFine(tMax(dNow - dLost, 0.0));
   m_UnconsolidatedSediment.SetNotchFineLost(0);

   dNow  = m_UnconsolidatedSediment.dGetSand();
   dLost = m_UnconsolidatedSediment.dGetNotchSandLost();
   m_UnconsolidatedSediment.SetSand(tMax(dNow - dLost, 0.0));
   m_UnconsolidatedSediment.SetNotchSandLost(0);

   dNow  = m_UnconsolidatedSediment.dGetCoarse();
   dLost = m_UnconsolidatedSediment.dGetNotchCoarseLost();
   m_UnconsolidatedSediment.SetCoarse(tMax(dNow - dLost, 0.0));
   m_UnconsolidatedSediment.SetNotchCoarseLost(0);

   dNow  = m_ConsolidatedSediment.dGetFine();
   dLost = m_ConsolidatedSediment.dGetNotchFineLost();
   m_ConsolidatedSediment.SetFine(tMax(dNow - dLost, 0.0));
   m_ConsolidatedSediment.SetNotchFineLost(0);

   dNow  = m_ConsolidatedSediment.dGetSand();
   dLost = m_ConsolidatedSediment.dGetNotchSandLost();
   m_ConsolidatedSediment.SetSand(tMax(dNow - dLost, 0.0));
   m_ConsolidatedSediment.SetNotchSandLost(0);

   dNow  = m_ConsolidatedSediment.dGetCoarse();
   dLost = m_ConsolidatedSediment.dGetNotchCoarseLost();
   m_ConsolidatedSediment.SetCoarse(tMax(dNow - dLost, 0.0));
   m_ConsolidatedSediment.SetNotchCoarseLost(0);
}

double CRWCellLayer::dGetFineUnconsolidatedThickness(void) const
{
   return m_UnconsolidatedSediment.dGetFine();
}

double CRWCellLayer::dGetFineConsolidatedThickness(void) const
{
   return m_ConsolidatedSediment.dGetFine();
}

double CRWCellLayer::dGetSandUnconsolidatedThickness(void) const
{
   return m_UnconsolidatedSediment.dGetSand();
}

double CRWCellLayer::dGetSandConsolidatedThickness(void) const
{
   return m_ConsolidatedSediment.dGetSand();
}

double CRWCellLayer::dGetCoarseUnconsolidatedThickness(void) const
{
   return m_UnconsolidatedSediment.dGetCoarse();
}

double CRWCellLayer::dGetCoarseConsolidatedThickness(void) const
{
   return m_ConsolidatedSediment.dGetCoarse();
}

double CRWCellLayer::dGetUnconsolidatedThickness(void) const
{
   return (m_UnconsolidatedSediment.dGetFine() + m_UnconsolidatedSediment.dGetSand() + m_UnconsolidatedSediment.dGetCoarse());
}

double CRWCellLayer::dGetConsolidatedThickness(void) const
{
   return (m_ConsolidatedSediment.dGetFine() + m_ConsolidatedSediment.dGetSand() + m_ConsolidatedSediment.dGetCoarse());
}

double CRWCellLayer::dGetTotalThickness(void) const
{
   return (m_UnconsolidatedSediment.dGetFine() + m_UnconsolidatedSediment.dGetSand() + m_UnconsolidatedSediment.dGetCoarse() + m_ConsolidatedSediment.dGetFine() + m_ConsolidatedSediment.dGetSand() + m_ConsolidatedSediment.dGetCoarse());
}

double CRWCellLayer::dGetNotchUnconsolidatedLost(void) const
{
   return (m_UnconsolidatedSediment.dGetNotchFineLost() + m_UnconsolidatedSediment.dGetNotchSandLost() + m_UnconsolidatedSediment.dGetNotchCoarseLost());
}

double CRWCellLayer::dGetNotchConsolidatedLost(void) const
{
   return (m_ConsolidatedSediment.dGetNotchFineLost() + m_ConsolidatedSediment.dGetNotchSandLost() + m_ConsolidatedSediment.dGetNotchCoarseLost());
}

// double CRWCellLayer::dGetVolSedFraction(void) const
// {
//    return m_VdolSedFraction;
// }

// void CRWCellLayer::SetVolSedFraction(double const dNewVolSedFraction)
// {
//    m_VdolSedFraction = dNewVolSedFraction;
// }
//
// double CRWCellLayer::dGetMechResistance(void) const
// {
//    return m_dMechResistance;
// }

// void CRWCellLayer::SetMechResistance(double const dNewMechResistance)
// {
//    m_dMechResistance = dNewMechResistance;
// }

// double CRWCellLayer::dGetConsolidationStatus(void) const
// {
//    return m_dConsolidationStatus;
// }

// void CRWCellLayer::SetConsolidationStatus(double const dNewConsolidationStatus)
// {
//    m_dConsolidationStatus = dNewConsolidationStatus;
// }

