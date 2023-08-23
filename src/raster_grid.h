/*!
 *
 * \class CGeomRasterGrid
 * \brief Geometry cass used to represent the raster grid of cell objects
 * \details TODO This is a more detailed description of the CGeomRasterGrid class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file raster_grid.h
 * \brief Contains CGeomRasterGrid definitions
 *
 */

#ifndef RASTERGRID_H
#define RASTERGRID_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "cell.h"


class CGeomCell;           // Forward declaration
class CSimulation;         // Ditto

class CGeomRasterGrid
{
   friend class CSimulation;
   friend class CGeomProfile;

private:
   double
      m_dD50Fine,
      m_dD50Sand,
      m_dD50Coarse;

   CSimulation* m_pSim;

   // NOTE seems to be faster than 2D STL vectors
   CGeomCell** m_Cell;

public:
   explicit CGeomRasterGrid(CSimulation*);
   ~CGeomRasterGrid(void);

   CSimulation* pGetSim(void);
//    CGeomCell* pGetCell(int const, int const);
   int nCreateGrid(void);
};
#endif // RASTERGRID_H
