#ifndef _PAIRTREE_H
#define _PAIRTREE_H

#include <fstream>
#include <cassert>

#include "eef1.hpp"

class CTerm
{
public:
  typedef enum {IDENTICAL, DIFFERENT} TYPE;

  CTerm(int depth, TYPE type);
  ~CTerm() 
  {
    for (int i = 0; i < 4; i++)
      if (m_pChildren[i])
	delete m_pChildren[i];
  }

  // Compute the electrostatic contribution of the two nodes
  static REAL computeElectrostatics(int type1, int type2,  
				    REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
				    int diff);

  // Compute the vdW energy between two links. 
  static REAL computeVdW(int type1, int type2, 
			 REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
			 int diff);

  // Compute the solvation energy between two links
  static REAL computeSolvation(int type1, int type2, 
			       REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE], 
			       int diff);

  // Recompute the energy stored at this node as the sum of the energy stored at
  // all its child nodes
  void recompute(bool b = true)
  {
    m_undoEnergy = m_energy;

    m_energy = m_pChildren[0]->m_energy +
      m_pChildren[1]->m_energy;

    if (b)
      {
	 m_energy += m_pChildren[2]->m_energy;
	if (m_pChildren[3])
	  m_energy += m_pChildren[3]->m_energy;
      }
  }

  // Set the energy to 0.
  void reset()
  { 
    m_undoEnergy = m_energy; 
    m_energy = 0.0; 
  }

  // Set the energy to a new value.
  void set(REAL energy)
  {
    m_undoEnergy = m_energy; 
    m_energy = energy; 
  }

  // Undo the change in energy caused by the latest move.
  void undo()
  {
    m_energy = m_undoEnergy; 
  }

  CTerm * getChild(int i)
  { return m_pChildren[i]; }

  CTerm * m_pChildren[4];
 
  REAL m_energy;
  REAL m_undoEnergy;
 };


#endif
