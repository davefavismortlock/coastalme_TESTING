/*!
 *
 * \file sediment_input_event.cpp
 * \brief CSedInputEvent routines
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

#include "cme.h"
#include "sediment_input_event.h"


//! Constructor with 8 parameters
CSedInputEvent::CSedInputEvent(int const nIDIn, unsigned long const ulTimeStepIn, double const dFineIn, double const dSandIn, double const dCoarseIn, double const dLenIn, double const dWidthIn): //, double const dThickIn):
   m_nLocationID(nIDIn),
   m_ulEventTimeStep(ulTimeStepIn),
   m_dFineSedVol(dFineIn),
   m_dSandSedVol(dSandIn),
   m_dCoarseSedVol(dCoarseIn),
   m_dLen(dLenIn),
   m_dWidth(dWidthIn)
   // m_dThick(dThickIn)
{
}

CSedInputEvent::~CSedInputEvent(void)
{
}

int CSedInputEvent::nGetLocationID(void)
{
   return m_nLocationID;
}

unsigned long CSedInputEvent::ulGetEventTimeStep(void)
{
   return m_ulEventTimeStep;
}

double CSedInputEvent::dGetFineSedVol(void)
{
   return m_dFineSedVol;
}

double CSedInputEvent::dGetSandSedVol(void)
{
   return m_dSandSedVol;
}

double CSedInputEvent::dGetCoarseSedVol(void)
{
   return m_dCoarseSedVol;
}

double CSedInputEvent::dGetLen(void)
{
   return m_dLen;
}

double CSedInputEvent::dGetWidth(void)
{
   return m_dWidth;
}

// double CSedInputEvent::dGetThick(void)
// {
//    return m_dThick;
// }

