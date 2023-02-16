/*!
 *
 * \file 2di_shape.cpp
 * \brief Abstract class, used as a base class for integer 2D objects (line, area, etc.)
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
#include <vector>

#include "2di_shape.h"

//! Constructor, no parameters
CA2DIShape::CA2DIShape(void)
{
}

//! Destructor
CA2DIShape::~CA2DIShape(void)
{
}


CGeom2DIPoint& CA2DIShape::operator[] (int const n)
{
   // NOTE No safety check
   return m_VPoints[n];
}


CGeom2DIPoint& CA2DIShape::Back(void)
{
   return m_VPoints.back();
}

vector<CGeom2DIPoint>* CA2DIShape::pPtiVGetPoints(void)
{
   return &m_VPoints;
}


void CA2DIShape::Clear(void)
{
   m_VPoints.clear();
}

void CA2DIShape::Resize(const int nSize)
{
   m_VPoints.resize(nSize);
}

int CA2DIShape::nGetSize(void) const
{
   return static_cast<int>(m_VPoints.size());
}


// void CA2DIShape::InsertAtFront(int const nX, int const nY)
// {
//    m_VPoints.insert(m_VPoints.begin(), CGeom2DIPoint(nX, nY));
// }

void CA2DIShape::Append(CGeom2DIPoint const* pPtiNew)
{
   m_VPoints.push_back(*pPtiNew);
}

void CA2DIShape::Append(int const nX, int const nY)
{
   m_VPoints.push_back(CGeom2DIPoint(nX, nY));
}

void CA2DIShape::AppendIfNotAlready(int const nX, int const nY)
{
   CGeom2DIPoint PtiIn(nX, nY);

   if (m_VPoints.empty())
      m_VPoints.push_back(PtiIn);

   else if (m_VPoints.back() != &PtiIn)
      m_VPoints.push_back(PtiIn);
}


// void CA2DIShape::SetPoints(const vector<CGeom2DIPoint>* VNewPoints)
// {
//    m_VPoints = *VNewPoints;
// }

// int CA2DIShape::nLookUp(CGeom2DIPoint* Pti)
// {
//    auto it = std::find(m_VPoints.begin(), m_VPoints.end(), *Pti);
//    if (it != m_VPoints.end())
//       return it - m_VPoints.begin();
//    else
//       return -1;
// }
