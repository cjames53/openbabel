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

using namespace std;

namespace OpenBabel
{

  OBNbrList::OBNbrList(OBMol* mol, double rcut, int boxSize) : m_mol(mol), m_rcut(rcut)
  {
    m_rcut2 = rcut*rcut;
    m_boxSize = boxSize;
    m_edgeLength = m_rcut / m_boxSize;
    m_updateCounter = 0;
    m_coords = m_mol->GetCoordinates();
 
    initOffsetMap();
    initOneTwo();
    initCells();
  }

  std::vector<OBAtom*> OBNbrList::nbrs(OBAtom *atom)
  {
    Eigen::Vector3d atomPos(m_coords + 3 * (atom->GetIdx()-1) );
    m_r2.clear();
    m_r2.reserve(m_mol->NumAtoms());
    std::vector<OBAtom*> atoms;
    atoms.reserve(m_mol->NumAtoms());
    Eigen::Vector3i index(cellIndexes(atomPos));

    std::vector<Eigen::Vector3i>::const_iterator i;
    // Use the offset map to find neighboring cells
    for (i = m_offsetMap.begin(); i != m_offsetMap.end(); ++i) {
      // add the offset to the cell index for atom's cell
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
        if (R2 > m_rcut2)
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
    m_oneTwo.resize(m_mol->NumAtoms());
    m_oneThree.resize(m_mol->NumAtoms());

    FOR_ATOMS_OF_MOL (atom, m_mol) {
      FOR_NBORS_OF_ATOM (nbr1, &*atom) {
        m_oneTwo[atom->GetIdx()-1].push_back(nbr1->GetIdx());
        m_oneTwo[nbr1->GetIdx()-1].push_back(atom->GetIdx());
            
        FOR_NBORS_OF_ATOM (nbr2, &*nbr1) {
          if (atom->GetIdx() == nbr2->GetIdx())
            continue;

          m_oneThree[atom->GetIdx()-1].push_back(nbr2->GetIdx());
          m_oneThree[nbr2->GetIdx()-1].push_back(atom->GetIdx());
        }
      }
    }
  }
      
  void OBNbrList::initCells()
  {
    // find min & max
    FOR_ATOMS_OF_MOL (atom, m_mol) {
      Eigen::Vector3d pos(m_coords + 3 * (atom->GetIdx()-1) );
  
      if (atom->GetIdx() == 1) {
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
    // the last cell is always empty and can be used for all ghost cells
    // in non-periodic boundary conditions.
    m_cells.resize(m_xyDim * m_dim.z() + 1);
    FOR_ATOMS_OF_MOL (atom, m_mol) {
      Eigen::Vector3d pos(m_coords + 3 * (atom->GetIdx()-1) );
      m_cells[cellIndex(pos)].push_back(&*atom);
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
