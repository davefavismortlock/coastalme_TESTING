/*!
 *
 * \class CSedInputEvent
 * \brief Class used to represent a sediment input event
 * \details This class represent a sediment input event such as sediment derived from inland (e.g. at the mouth of a gully or rambla), or sediment from an intervention such as beach nourishment
 * \author David Favis-Mortlock
 * \author Andres Payo
 * \date 2023
 * \copyright GNU General Public License
 *
 * \file sediment_input_event.h
 * \brief Contains CSedInputEvent definitions
 *
 */

#ifndef CSEDINPUT_H
#define CSEDINPUT_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
class CSedInputEvent
{
private:
   int m_nLocationID;        // The location ID in the shapefile
   unsigned long m_ulEventTimeStep;        // The timing of the sediment input event
   double
      m_dFineSedVol,         // Volume (m3) of fine sediment
      m_dSandSedVol,         // Volume (m3) of sand-sized sediment
      m_dCoarseSedVol,       // Volume (m3) of coarse sediment
      m_dLen,                // The coast-normal length (m) of the sediment block
      m_dWidth;              // The along-coast width (m) of the sediment block
      // m_dThick;              // The along-coast thickness (m) of the sediment block

public:
   CSedInputEvent(int const, unsigned long const, double const, double const, double const, double const, double const); //, double const);
   ~CSedInputEvent(void);

   int nGetLocationID(void);
   unsigned long ulGetEventTimeStep(void);
   double dGetFineSedVol(void);
   double dGetSandSedVol(void);
   double dGetCoarseSedVol(void);
   double dGetLen(void);
   double dGetWidth(void);
   // double dGetThick(void);
};
#endif // CSEDINPUT_H
