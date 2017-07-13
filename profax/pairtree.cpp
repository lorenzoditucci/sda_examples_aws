#include <fstream>
#include <cassert>

#include "pairtree.hpp"
#include "eef1.hpp"

// Create a node of the energytree. A recursive procedure that
// creates all of the energytree.
CTerm::CTerm(int depth, TYPE type) : 
  m_energy(0.0), m_undoEnergy(0.0)
{
  if (depth > 0)
    {
      if (type == IDENTICAL)
	{
	  m_pChildren[0] = new CTerm(depth-1, IDENTICAL);
	  m_pChildren[1] = new CTerm(depth-1, DIFFERENT);
	  m_pChildren[2] = new CTerm(depth-1, IDENTICAL);
	  m_pChildren[3] = NULL;
	} 
      else // TYPE == DIFFERNT
	{
	  m_pChildren[0] = new CTerm(depth-1, DIFFERENT);
	  m_pChildren[1] = new CTerm(depth-1, DIFFERENT);
	  m_pChildren[2] = new CTerm(depth-1, DIFFERENT);
	  m_pChildren[3] = new CTerm(depth-1, DIFFERENT);
	}
    }
  else
    {
      m_pChildren[0] = m_pChildren[1] = m_pChildren[2] = 
	m_pChildren[3] = NULL;
    }
}


// Compute the electrostatic contribution of two leaf nodes
// from the matrix of distances between the atoms.
REAL CTerm::computeElectrostatics(int type1, int type2,  
				  REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
				  int diff)
{
  assert(diff >= 0);
  
  int size1 = SIDECHAIN::m_aalist[type1]->m_size;
  int size2 = SIDECHAIN::m_aalist[type2]->m_size;
  
  // Proline is treated differently because the Cd atom is in a group
  // with the Ca atom which is in a different leaf node.
  if (type1 == BBP)
    size1++;
  else if (type1 == PRO)
    size1--;
  
  if (type2 == BBP)
    size2++;
  else if (type2 == PRO)
    size2--;
  
  int ngr1 = SIDECHAIN::m_aalist[type1]->m_nGroups;
  int ngr2 = SIDECHAIN::m_aalist[type2]->m_nGroups;
  
  int * groups1 = SIDECHAIN::m_aalist[type1]->m_groups;
  int * groups2 = SIDECHAIN::m_aalist[type2]->m_groups;
  
  REAL * charges1 = SIDECHAIN::m_aalist[type1]->m_charges;
  REAL * charges2 = SIDECHAIN::m_aalist[type2]->m_charges;
  
  REAL sum = 0.0;

  for (int g1 = 0; g1 < ngr1; g1++)
    for (int g2 = 0; g2 < ngr2; g2++)
      {
	int end1, end2, start1, start2;
	
	start1 = groups1[g1];
	start2 = groups2[g2];
	
	if (g1 == ngr1 - 1)
	  end1 = size1;
	  else 
	    end1 = groups1[g1 + 1];
	
	if (g2 == ngr2 - 1)
	  end2 = size2;
	else 
	  end2 = groups2[g2 + 1];
	
	// Decide whether the two groups are close enough to interact.
	// They are close enough if they contain an interacting pair 
	// of atoms.
	bool bCompute = false;
	for (int i = start1; i < end1; i++)
	  for(int j = start2; j < end2; j++)
	    {
	      if (dists[i][j] < CUTOFF_DISTANCE_2)
		{
		  bCompute = true;
		  i = 1000;
		  break;
		}
	    }
	
	// If the groups are interacting, compute their contribution.
	if (bCompute)
	  {
	    for (int i = start1; i < end1; i++)
	      for(int j = start2; j < end2; j++)
		{
		  int ex = isExcluded(diff, type1, i, type2, j);
		  if (ex != EXCLUDED)
		  {
		    REAL w = compute_ES(charges1[i], charges2[j],
					dists[i][j]);

		    if (ex == PAIR1_4)
		      w *= 0.4;

		    sum += w;
		  }
		}
	  }
      }
  
  return sum;
}

// Compute the vdW energy between two leaf nodes. 
REAL CTerm::computeVdW(int type1, int type2, 
		       REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE],
		       int diff)
{
  assert(diff >= 0);
  
  int size1 = SIDECHAIN::m_aalist[type1]->m_size;
  int size2 = SIDECHAIN::m_aalist[type2]->m_size;
  
  REAL sum = 0.0;
  
  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      {
	int ex = isExcluded(diff, type1, i, type2, j);
	if (ex != EXCLUDED)
	  {
	    sum += compute_vdW(SIDECHAIN::m_aalist[type1]->m_aTypes[i],
			       SIDECHAIN::m_aalist[type2]->m_aTypes[j],
			       dists[i][j], ex == PAIR1_4);
	  }
      }
  
  return sum;
}

// Compute the solvation energy between two links
REAL CTerm::computeSolvation(int type1, int type2, 
			     REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE], 
			     int diff)
{
  assert(diff >= 0);
  
  int size1 = SIDECHAIN::m_aalist[type1]->m_size;
  int size2 = SIDECHAIN::m_aalist[type2]->m_size;
  
  REAL sum = 0.0;
  
  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      {
	int ex = isExcluded(diff, type1, i, type2, j);
	if (ex != EXCLUDED)
	  {
	    REAL lambda1 = getLambda(type1, i, SIDECHAIN::m_aalist[type1]->m_aTypes[i]);
	    REAL lambda2 = getLambda(type2, j, SIDECHAIN::m_aalist[type2]->m_aTypes[j]);
	    sum += computeSolventEffect(SIDECHAIN::m_aalist[type1]->m_aTypes[i],
					SIDECHAIN::m_aalist[type2]->m_aTypes[j],
					dists[i][j], lambda1, lambda2, ex == PAIR1_4);
	  }	
      }
  
  return sum;
}

