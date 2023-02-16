/*!
 *
 * \class CRWCellSediment
 * \brief Real-world class used to represent the sediment (either consolidated or unconsolidated) associated with a cell layer object
 * \details TODO This is a more detailed description of the CRWCellSediment class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file cell_sediment.h
 * \brief Contains CRWCellSediment definitions
 *
 */

#ifndef SEDIMENT_H
#define SEDIMENT_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
class CRWCellSediment
{
private:
   double
      m_dFine,                // Depth equivalent in m
      m_dNotchFineLost,       // Depth equivalent (m) of sediment lost via notch incision
      m_dSand,                // Depth equivalent in m
      m_dNotchSandLost,       // Depth equivalent (m) of sediment lost via notch incision
      m_dCoarse,              // Depth equivalent in m
      m_dNotchCoarseLost,     // Depth equivalent (m) of sediment lost via notch incision
      m_dTotSedimentInput;    // Total depth (m, from all size classes) of sediment input events, since start of simulation

public:
   CRWCellSediment(void);
   CRWCellSediment(CRWCellSediment const&);           // Copy constructor defined explicitly, to stop cppcheck from complaining

   CRWCellSediment& operator= (const CRWCellSediment&);

   void SetFine(double const);
   double dGetFine(void) const;
   void AddFine(double const);

   void SetSand(double const);
   double dGetSand(void) const;
   void AddSand(double const);

   void SetCoarse(double const);
   double dGetCoarse(void) const;
   void AddCoarse(double const);

   void SetNotchFineLost(double const);
   void IncrNotchFineLost(double const);
   double dGetNotchFineLost(void) const;

   void SetNotchSandLost(double const);
   void IncrNotchSandLost(double const);
   double dGetNotchSandLost(void) const;

   void SetNotchCoarseLost(double const);
   void IncrNotchCoarseLost(double const);
   double dGetNotchCoarseLost(void) const;

   void AddToTotSedimentInputDepth(double const);
   double dGetTotSedimentInputDepth(void) const;
};
#endif // SEDIMENT_H
