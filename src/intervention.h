/*!
 *
 * \class CRWIntervention
 * \brief Real-world class used to represent the 'intervention' category of coastal landform objects
 * \details TODO This is a more detailed description of the CRWIntervention class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2021
 * \copyright GNU General Public License
 *
 * \file intervention.h
 * \brief Contains CRWIntervention definitions
 *
 */

#ifndef INTERVENTION_H
#define INTERVENTION_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "coast.h"
#include "coast_landform.h"


class CRWIntervention : public CACoastLandform
{
private:

public:
   CRWIntervention(CRWCoast*, int const, int const);
   ~CRWIntervention(void);

   void Display(void) override;           // Instantiates the pure virtual function in the abstract parent class, so that CRWIntervention is not an abstract class. But this instatiation never gets called, which seems like a waste of time. Alternative?
};
#endif // INTERVENTION_H

