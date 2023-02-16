/*!
 *
 * \brief Definitions of some routines from the linear interp library
 * \details TODO This is a more detailed description of the linear interp routines.
 * \author Modified by Andres Payo and David Favis-Mortlock
 * \date 2023
 * \copyright GNU Lesser General Public License
 *
 * \file linearinterp.cpp
 * \brief Contains definitions of routines from the linear interp library
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include "linearinterp.h"

int nNearestNeighbourIndex (vector<double> const* pVdX, double const dValue)
{
   double dDist = DBL_MAX;
   int nIdx = 0;

   for (unsigned int i = 0; i < pVdX->size(); ++i)
   {
      double dNewDist = tAbs (dValue - pVdX->at (i));

      if (dNewDist <= dDist)
      {
         dDist = dNewDist;
         nIdx = i;
      }
   }

   return nIdx;
}


vector<double> VdInterp1 (vector<double> const* pVdX, vector<double> const* pVdY, vector<double> const* pVdX_new)
{
   vector<double> VdY_new;
   double dX, dY;
   unsigned int
   x_max_idx = static_cast<unsigned int> (pVdX->size()) - 1,
   x_new_size = static_cast<unsigned int> (pVdX_new->size());

   VdY_new.reserve (x_new_size);

   for (unsigned int i = 0; i < x_new_size; ++i)
   {
      unsigned int idx = nNearestNeighbourIndex (pVdX, pVdX_new->at (i));

      if (pVdX->at (idx) > pVdX_new->at (i))
      {
         if (idx > 0)
         {
            dX = pVdX->at (idx) - pVdX->at (idx - 1);
            dY = pVdY->at (idx) - pVdY->at (idx - 1);
         }
         else
         {
            dX = pVdX->at (idx + 1) - pVdX->at (idx);
            dY = pVdY->at (idx + 1) - pVdY->at (idx);
         }
      }
      else
      {
         if (idx < x_max_idx)
         {
            dX = pVdX->at (idx + 1) - pVdX->at (idx);
            dY = pVdY->at (idx + 1) - pVdY->at (idx);
         }
         else
         {
            dX = pVdX->at (idx) - pVdX->at (idx - 1);
            dY = pVdY->at (idx) - pVdY->at (idx - 1);
         }
      }

      double m = dY / dX;
      double b = pVdY->at (idx) - pVdX->at (idx) * m;

      VdY_new.push_back (pVdX_new->at (i) * m + b);
   }

   return VdY_new;
}
