/*********************************************************************
nbrlist.cpp - NbrList class

Copyright (C) 2009 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "obnbrlist.h"

namespace OpenBabel
{

  OBNbrList::OBNbrList(const std::vector<OBAtom*> atoms, double rcut, int boxSize) : m_atoms(atoms), m_rcut(rcut)
  {
    m_boxSize = boxSize;
    m_edgeLength = m_rcut / m_boxSize;
    m_updateCounter = 0;
 
    initOffsetMap();
    initOneTwo();
    initCells();
  }

  std::vector<OBAtom*> OBNbrList::nbrs(OBAtom *atom)
  {
    Eigen::Vector3d atomPos(atom->GetVector().AsArray());
    m_r2.clear();
    m_r2.reserve(m_atoms.size() / m_xyDim);
    std::vector<OBAtom*> atoms;
    atoms.reserve(m_atoms.size() / m_xyDim);
    Eigen::Vector3i index(cellIndexes(Eigen::Vector3d(atom->GetVector().AsArray())));
    const double rcut2 = m_rcut * m_rcut;

    std::vector<Eigen::Vector3i>::const_iterator i;
    for (i = m_offsetMap.begin(); i != m_offsetMap.end(); ++i) {
      Eigen::Vector3i offset = index + *i;
 
      if (offset.x() < 0) continue;
      if (offset.y() < 0) continue;
      if (offset.z() < 0) continue;
      if (offset.x() >= m_dim.x()) continue;
      if (offset.y() >= m_dim.y()) continue;
      if (offset.z() >= m_dim.z()) continue;

      unsigned int cell = cellIndex(offset);
      for (atom_iter j = m_cells[cell].begin(); j != m_cells[cell].end(); ++j) {
        // make sure to only return unique pairs
        if (atom->GetIdx() >= (*j)->GetIdx())
          continue;
        if (IsOneTwo(atom->GetIdx(), (*j)->GetIdx()))
          continue;
        if (IsOneThree(atom->GetIdx(), (*j)->GetIdx()))
          continue;

        const double R2 =  (Eigen::Vector3d((*j)->GetVector().AsArray())-atomPos).squaredNorm();
        if (R2 > rcut2)
          continue;

        m_r2.push_back(R2);
        atoms.push_back(*j);
      }
    }

    return atoms;
  }

  void OBNbrList::update()
  {
    m_updateCounter++;

    if (m_updateCounter > 10) {
      initCells();    
      m_updateCounter = 0;
    }
  }

  void OBNbrList::initOneTwo()
  {
    m_oneTwo.resize(m_atoms.size());
    m_oneThree.resize(m_atoms.size());

    for (atom_iter a = m_atoms.begin(); a != m_atoms.end(); ++a) {
      FOR_NBORS_OF_ATOM (nbr1, *a) {
        m_oneTwo[(*a)->GetIdx()-1].push_back(nbr1->GetIdx());
        m_oneTwo[nbr1->GetIdx()-1].push_back((*a)->GetIdx());
            
        FOR_NBORS_OF_ATOM (nbr2, &*nbr1) {
          if ((*a)->GetIdx() == nbr2->GetIdx())
            continue;

          m_oneThree[(*a)->GetIdx()-1].push_back(nbr2->GetIdx());
          m_oneThree[nbr2->GetIdx()-1].push_back((*a)->GetIdx());
        }
      }
    }
  }
      
  void OBNbrList::initCells()
  {
    // find min & max
    for (atom_iter i = m_atoms.begin(); i != m_atoms.end(); ++i) {
      Eigen::Vector3d pos((*i)->GetVector().AsArray());
  
      if ((*i)->GetIdx() == 1) {
        m_min = m_max = pos;
      } else {
        if (pos.x() > m_max.x())
          m_max.x() = pos.x();
        else if (pos.x() < m_min.x())
          m_min.x() = pos.x();

        if (pos.y() > m_max.y())
          m_max.y() = pos.y();
        else if (pos.y() < m_min.y())
          m_min.y() = pos.y();
 
        if (pos.z() > m_max.z())
          m_max.z() = pos.z();
        else if (pos.z() < m_min.z())
          m_min.z() = pos.z();
      }
    }
 
    // set the dimentions
    m_dim.x() = floor( (m_max.x() - m_min.x()) /  m_edgeLength) + 1;
    m_dim.y() = floor( (m_max.y() - m_min.y()) /  m_edgeLength) + 1;
    m_dim.z() = floor( (m_max.z() - m_min.z()) /  m_edgeLength) + 1;
    m_xyDim = m_dim.x() * m_dim.y();
      
    // add atoms to their cells
    m_cells.clear();
    m_cells.resize(m_xyDim * m_dim.z());
    for (atom_iter i = m_atoms.begin(); i != m_atoms.end(); ++i) {
      Eigen::Vector3d pos((*i)->GetVector().AsArray());
      m_cells[cellIndex(pos)].push_back(*i);
    }
      
  }

  void OBNbrList::initOffsetMap()
  {
    int dim = 2 * m_boxSize + 1;
    m_offsetMap.resize(dim * dim * dim);
    for (int i = 0; i < dim; ++i) 
      for (int j = 0; j < dim; ++j) 
        for (int k = 0; k < dim; ++k) 
          m_offsetMap[offsetIndex(i,j,k)] = Eigen::Vector3i(i - m_boxSize, j - m_boxSize, k - m_boxSize);
  }

} // end namespace OpenBabel

//! \file nbrlist.cpp
//! \brief NbrList class
