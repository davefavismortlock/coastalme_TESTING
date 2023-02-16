/*!
 *
 * \file locate_estuaries.cpp
 * \brief Finds estuaries on the raster grid
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include "cme.h"
#include "simulation.h"


/*===============================================================================================================================

 Locate estuaries

===============================================================================================================================*/
int CSimulation::nLocateAllEstuaries(void)
{
   // TODO
   // For each estuary (if any):
   //   FindChannelNetwork();
   //   GetEstuaryEvents(Area, Volume);

   //   For each channel section:
   //      CreateEstuaryCrossSection();
   //      GetSectionsEvents(Geometry, Properties);

   return RTN_OK;
}
