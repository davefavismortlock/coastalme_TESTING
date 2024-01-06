/*!
 *
 * \file interpolate.cpp
 * \brief Linear interpolation routines
 * \details TODO A more detailed description of these routines.
 * \author http://www.cplusplus.com/forum/general/216928/
 * \author Modified by David Favis-Mortlock and Andres Payo
 * \date 2024
 * \copyright GNU Lesser General Public License
 *
 */

/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "cme.h"

//===============================================================================================================================
//! Returns interpolated value at x from parallel arrays (VdXdata, VdYdata). Assumes that VdXdata has at least two elements, is sorted and is strictly monotonically increasing. The boolean argument extrapolate determines behaviour beyond ends of array (if needed). For this version, both lots of data are doubles
//===============================================================================================================================
double dInterpolate(vector<double> const* pVdXdata, vector<double> const* pVdYdata, double dX, bool bExtrapolate)
{
   int size = static_cast<int>(pVdXdata->size());

   int i = 0;                                   // Find left end of interval for interpolation
   if (dX >= pVdXdata->at(size - 2))            // Special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while (dX > pVdXdata->at(i+1)) i++;
   }
   
   double 
      dXL = pVdXdata->at(i),
      dYL = pVdYdata->at(i),
      dXR = pVdXdata->at(i+1),
      dYR = pVdYdata->at(i+1);                  // Points on either side (unless beyond ends)
      
   if (! bExtrapolate)                          // If beyond ends of array and not extrapolating
   {
      if (dX < dXL) dYR = dYL;
      if (dX > dXR) dYL = dYR;
   }

   double ddYdX = (dYR - dYL) / (dXR - dXL);    // Gradient

   return (dYL + ddYdX * (dX - dXL));           // Linear interpolation
}

//===============================================================================================================================
//! Returns interpolated value at x from parallel arrays (VdXdata, VdYdata). Assumes that VdXdata has at least two elements, is sorted and is strictly monotonically increasing. The boolean argument extrapolate determines behaviour beyond ends of array (if needed). For this version, one lot of data is integer and the other is double
//===============================================================================================================================
double dInterpolate(vector<int> const* pVnXdata, vector<double> const* pVdYdata, int nX, bool bExtrapolate )
{
   unsigned int nSize = static_cast<unsigned int>(pVnXdata->size());

   int i = 0;                                   // Find left end of interval for interpolation
   if (nX >= pVnXdata->at(nSize - 2))           // Special case: beyond right end
   {
      i = nSize - 2;
   }
   else
   {
      while (nX > pVnXdata->at(i+1)) i++;
   }
   
   int 
      nXL = pVnXdata->at(i),
      nXR = pVnXdata->at(i+1);
      
   double 
      dYL = pVdYdata->at(i),
      dYR = pVdYdata->at(i+1);                  // Points on either side (unless beyond ends)
      
   if (! bExtrapolate)                          // If beyond ends of array and not extrapolating
   {
      if (nX < nXL) dYR = dYL;
      if (nX > nXR) dYL = dYR;
   }

   double ddYdX = (dYR - dYL) / static_cast<double>(nXR - nXL);      // Gradient

   return dYL + ddYdX * static_cast<double>(nX - nXL);               // Linear interpolation
}

//======================================================================
//int main()
//{
//   // Original data
//   vector<double> VdXdata = { 1, 5, 10, 15, 20 };
//   vector<double> VdYdata = { 0.3, 0.5, 0.8, 0.1, 0.14 };

//   // Set up some points for interpolation in xVals
//   const int NPTS = 20;
//   vector<double> xVals, yVals;
//   for ( int i = 1; i <= NPTS; i++ ) xVals.push_back( (double)i );

//   // Interpolate
//   for ( double x : xVals )
//   {
//      double y = dInterpolate( VdXdata, VdYdata, x, true );
//      yVals.push_back( y );
//   }
//
//   // Output
//   #define SP << fixed << setw( 15 ) << setprecision( 6 ) <<
//   #define NL << '\n'
//   cout << "Original data:\n";
//   for ( int i = 0; i < VdXdata.size(); i++ ) cout SP VdXdata[i] SP VdYdata[i] NL;
//   cout << "\nInterpolated data:\n";
//   for ( int i = 0; i < xVals.size(); i++ ) cout SP xVals[i] SP yVals[i] NL;
//}

