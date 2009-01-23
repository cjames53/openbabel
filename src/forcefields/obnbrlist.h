/*********************************************************************
nbrlist.h - NbrList class

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

/**
 * @class OBNbrList
 *
 * Based on:
 * Near-neighbor calculations using a modified cell-linked list method
 * Mattson, W.; B. M. Rice (1999). "Near-neighbor calculations using a 
 * modified cell-linked list method". Computer Physics Communications 
 * 119: 135. 
 *
 * http://dx.doi.org/10.1016/S0010-4655%2898%2900203-3
 *
 */

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include <Eigen/Core>

namespace OpenBabel
{

  class OBNbrList
  {
    private:
      typedef std::vector<OBAtom*>::iterator atom_iter;

    public:
      OBNbrList(OBMol *mol, double rcut, int boxSize = 1);
      void update();
      
      std::vector<OBAtom*> nbrs(OBAtom *atom);
      inline double r2(unsigned int index) const
      {
        return m_r2.at(index);
      }

    private:
      inline unsigned int offsetIndex(int i, int j, int k) const
      {
        int dim = 2 * m_boxSize + 1;
        return i + j * dim + k * dim * dim;
      }

      inline unsigned int cellIndex(int i, int j, int k) const
      {
        return i + j * m_dim.x() + k * m_xyDim;
      }

      inline unsigned int cellIndex(const Eigen::Vector3i &index) const
      {
        return index.x() + index.y() * m_dim.x() + index.z() * m_xyDim;
      }


      inline unsigned int cellIndex(const Eigen::Vector3d &pos) const
      {
        return cellIndex( floor( (pos.x() - m_min.x()) / m_edgeLength ),
                          floor( (pos.y() - m_min.y()) / m_edgeLength ),
                          floor( (pos.z() - m_min.z()) / m_edgeLength ) );
      }

      inline Eigen::Vector3i cellIndexes(const Eigen::Vector3d &pos) const
      {
        Eigen::Vector3i index;
        index.x() = floor( (pos.x() - m_min.x()) / m_edgeLength );
        index.y() = floor( (pos.y() - m_min.y()) / m_edgeLength );
        index.z() = floor( (pos.z() - m_min.z()) / m_edgeLength );
        return index;
      }

      inline bool IsOneTwo(unsigned int i, unsigned int j) const 
      {
        --i;
        std::vector<unsigned int>::const_iterator iter;
        for (iter = m_oneTwo.at(i).begin(); iter != m_oneTwo.at(i).end(); ++iter)
          if (*iter == j)
            return true;

        return false;
      }
 
      inline bool IsOneThree(unsigned int i, unsigned int j) const
      {
        --i;
        std::vector<unsigned int>::const_iterator iter;
        for (iter = m_oneThree.at(i).begin(); iter != m_oneThree.at(i).end(); ++iter)
          if (*iter == j)
            return true;

        return false;
      }

      void initOneTwo();
      void initCells();
      void initOffsetMap();

      OBMol                              *m_mol;
      double                             *m_coords;
      std::vector<Eigen::Vector3i>        m_offsetMap;
      double                              m_rcut, m_rcut2;
      double                              m_edgeLength;
      int                                 m_boxSize;
      int                                 m_updateCounter;
    
      Eigen::Vector3d                     m_min, m_max;
      Eigen::Vector3i                     m_dim;
      double                              m_xyDim;
      std::vector<std::vector<OBAtom*> >  m_cells;

      std::vector<double>                 m_r2;
      
      std::vector<std::vector<unsigned int> > m_oneTwo;
      std::vector<std::vector<unsigned int> > m_oneThree;
  };
 
} // end namespace OpenBabel

//! \file nbrlist.h
//! \brief NbrList class
