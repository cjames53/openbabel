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

  /**
   *  
   *
   *
   *
   */
  class OBNbrList
  {
    private:
      typedef std::vector<OBAtom*>::iterator atom_iter;

    public:
      /**
       * Constructor.
       * @param mol The molecule containing the atoms
       * @param rcut The cut-off distance.
       * @param boxSize The number of cells per rcut distance.
       */
      OBNbrList(OBMol *mol, double rcut, int boxSize = 1);
      /**
       * Update the cells. While minimizing or running MD simulations,
       * atoms move and can go from on cell into the next. This function
       * should be called every 10-20 iterations to make sure the cells
       * stay accurate.
       */
      void update();
      /**
       * Get the near-neighbor atoms for @p atom. The squared distance is
       * checked and is cached for later use (see r2() function). 
       *
       * Note: Atoms in relative 1-2 and 1-3 positions are not returned.
       * The @p atom itself isn't added to the list.
       *
       * @param atom The atom for which to return the near-neighbors
       * @return The near-neighbors for @p atom
       */
      std::vector<OBAtom*> nbrs(OBAtom *atom);
      /**
       * Get the cached squared distance from the atom last used to call 
       * nbrs to the atom with @p index in the returned vector.
       * @param index The index for the atom in the vector of atoms returned by nbrs().
       * @return The cached squared distance.
       */
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

      inline unsigned int ghostIndex(int i, int j, int k) const
      {
        i += m_boxSize + 1;
        j += m_boxSize + 1;
        k += m_boxSize + 1;
        int xDim = 2 * m_boxSize + m_dim.x() + 2;
        int yDim = 2 * m_boxSize + m_dim.y() + 2;
        return i + j * xDim + k * xDim * yDim;
      }

      inline unsigned int ghostIndex(const Eigen::Vector3i &index) const
      {
        return ghostIndex(index.x(), index.y(), index.z());
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

      /**
       * @param i Index for the first atom (indexed from 1)
       * @param j Index for the first atom (indexed from 1)
       * @return True if atoms with index @p i and @p j are bonded (1-2)
       */
      inline bool IsOneTwo(unsigned int i, unsigned int j) const 
      {
        --i;
        std::vector<unsigned int>::const_iterator iter;
        for (iter = m_oneTwo.at(i).begin(); iter != m_oneTwo.at(i).end(); ++iter)
          if (*iter == j)
            return true;

        return false;
      }
 
      /**
       * @param i Index for the first atom (indexed from 1)
       * @param j Index for the first atom (indexed from 1)
       * @return True if atoms with index @p i and @p j are in a 1-3 position
       */
      inline bool IsOneThree(unsigned int i, unsigned int j) const
      {
        --i;
        std::vector<unsigned int>::const_iterator iter;
        for (iter = m_oneThree.at(i).begin(); iter != m_oneThree.at(i).end(); ++iter)
          if (*iter == j)
            return true;

        return false;
      }

      /**
       * Initialize the 1-2 and 1-3 cache.
       */
      void initOneTwo();
      void initCells();
      void initOffsetMap();
      void initGhostMap(bool periodic = false);
      bool insideShpere(const Eigen::Vector3i &index);

      OBMol                              *m_mol;
      double                             *m_coords;
      std::vector<Eigen::Vector3i>        m_offsetMap;
      std::vector<Eigen::Vector3i>        m_ghostMap;
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
