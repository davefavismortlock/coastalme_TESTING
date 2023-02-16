/*!
 *
 * \file utils_global.cpp
 * \brief Globally-available utility routines
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <cmath>
#include <cfloat>

#include <sstream>
using std::stringstream;

#include <iomanip>
using std::setw;

#include "cme.h"


/*==============================================================================================================================

 Correctly rounds doubles

==============================================================================================================================*/
double dRound(double const d)
{
   // Rounds positive or negative doubles correctly
   return ((d < 0.0) ? ceil(d - 0.5) : floor(d + 0.5));
}


/*==============================================================================================================================

 Version of the above that returns an int

 ==============================================================================================================================*/
int nRound(double const d)
{
   // Rounds positive or negative doubles correctly
   return static_cast<int>((d < 0.0) ? ceil(d - 0.5) : floor(d + 0.5));
}


// bool bIsWhole(double d)
// {
//    // From http://answers.yahoo.com/question/index?qid=20110320132617AAMdb7u
//    return (static_cast<int>(d) == d);
// }


/*==============================================================================================================================

 Checks to see if a string can be read as a valid double number. Does not find trailing (i.e.post-number) rubbish, but then neither does strtod(). From https://stackoverflow.com/questions/392981/how-can-i-convert-string-to-double-in-c

==============================================================================================================================*/
bool bIsStringValidDouble(string& str)
{
   std::istringstream iStr(str);
   double dDummy;

   if (!(iStr >> dDummy))
      return false;

   return true;
}


/*==============================================================================================================================

 Checks to see if a string can be read as a valid integer, from https://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int

==============================================================================================================================*/
bool bIsStringValidInt(string& str)
{
   // Trim leading whitespace
   size_t nPos = str.find_first_not_of(" \t");
   if (nPos != string::npos)
      str = str.substr(nPos);

   // If the first character is the sign, remove it
   if ((str[0] == '-') || (str[0] == '+'))
      str.erase(0, 1);

   // Now check that the string contains only numbers
   return (str.find_first_not_of("0123456789") == string::npos);
}


/*==============================================================================================================================

 Operator that inserts a given fill character, to a given width, into an output stream. From http://stackoverflow.com/questions/2839592/equivalent-of-02d-with-stdstringstream

==============================================================================================================================*/
ostream& operator<< (ostream& ostr, const FillToWidth& args)
{
   ostr.fill(args.chFill);
   ostr.width(args.nWidth);

   return ostr;
}

