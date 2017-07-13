#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <cfloat>
#include <fstream>
#include <iostream>

#include "pdb.hpp"
#include "chain.hpp"
#include "leaf.hpp"
#include "pairtree.hpp"
#include "eef1.hpp"



// Create a chain from a file describing the protein's structure.
CChain::CChain(const char * fname, int type)
{
  vector<int> AAtypes;
  vector<REAL> phis, psis;
  vector<int> rots;
  
  if (type == ANGS_FILE)
    load_angs(fname, AAtypes, phis, psis, rots);
  else
    {
      cout << "Bad type of input file!" << endl;
      exit(0);
    }
  create(AAtypes, phis, psis, rots);

  // Compute the reference solvent energy and the torsion energy for this chain.
  computeTorsionEnergy();
  m_solventE = computeSolventRefEnergy();
}

// Load the chain from a description that specifies phi/psi angles and rotamer indices.
void CChain::load_angs(const char * fname, vector<int> & AAtypes, vector<REAL> & phis,
		       vector<REAL> & psis, vector<int> & rots)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open input file " << fname << endl;
      exit(0);
    }

  char buf[100];

  float phi, psi;
  int rotIndex;
  char aa[5];
  while (!fin.eof())
    {
      fin.getline(buf, 99);

      if (strlen(buf) == 0)
	continue;

      int res = sscanf(buf, "%s %f %f %d", aa, &phi, &psi, &rotIndex);
      if (res < 1)
	{
	  cout << "File formatting error for .angs file: " << buf << endl;
	  exit(0);
	}

      AAtypes.push_back(getAA(aa));

      if (res >= 2)
	phis.push_back((REAL) (phi*M_PI/180));
      else
	phis.push_back(0.0);

      if (res >= 3)
	psis.push_back((REAL) (psi*M_PI/180));
      else
	psis.push_back(0.0);

      if (res == 4)
	{
	  rots.push_back(rotIndex);
	  if (rotIndex >= SIDECHAIN::m_aalist[AAtypes[AAtypes.size()-1]]->m_nRotamers)
	    {
	      cout << "Rotamer index " << rotIndex << " is too big for " << aa 
		   << rots.size() << " which has only " 
		   << SIDECHAIN::m_aalist[AAtypes[AAtypes.size()-1]]->m_nRotamers 
		   << " rotamers (the index is zero-based)!!!" << endl;
	      exit(0);
	    }
	}
      else
	rots.push_back(0);

      int i = rots.size()-1;
    }

  // The first phi and last psi are always 0
  phis[0] = 0.0; psis[psis.size()-1] = 0.0;

  cout << "loaded " << AAtypes.size() << " AAs" << endl;
  fin.close();
}


// Create the kinematic chain representaion of the protein from internal 
// coordinates (torsion angles).
void CChain::create(const vector<int> AAtypes, 
		    const vector<REAL> & phis, const vector<REAL> & psis, 
		    const vector<int> & rots)
{
  assert(AAtypes.size() > 1);

  const ROTAMER * rot;
  CLeaf * pLeaf1, * pLeaf2;

  pLeaf1 = new CLeaf(ROTAMER::UP, 0, &(D_PHI_AXIS[0]), NTR, 0);
  VcV(pLeaf1->m_translate, N_C_NTR_TRANS); 
  pLeaf1->setPrev(NULL);
  pLeaf1->rotate(phis[0]);
  m_links.push_back(pLeaf1);

  CLeaf * last = pLeaf1;
 
  for (int i = 0; i < AAtypes.size(); i++)
    {
      int bbn_type;
      if (i == AAtypes.size() - 1)
	bbn_type = CTR;
      else if (AAtypes[i+1] == PRO)
	bbn_type = BBP;
      else if (AAtypes[i+1] == GLY)
	bbn_type = BBG;
      else
	bbn_type = BBN;

      if (i%2 == 0)
	{
	  pLeaf1 = new CLeaf(ROTAMER::DOWN, rots[i], D_PSI_AXIS, AAtypes[i], 2*i+1);
	  VcV(pLeaf1->m_translate, D_CA_C_TRANS);
	  pLeaf1->rotate(psis[i]);

	  pLeaf2 = new CLeaf(ROTAMER::DOWN, 0, U_PHI_AXIS, bbn_type, 2*i+2);
	  VcV(pLeaf2->m_translate, UD_C_CA_TRANS);
	  if (bbn_type == BBP)
	    pLeaf2->rotate(-M_PI/3.0);
	  else
	    pLeaf2->rotate(phis[i+1]);
	}
      else
	{
	  pLeaf1 = new CLeaf(ROTAMER::UP, rots[i], U_PSI_AXIS, AAtypes[i], 2*i+1);
	  VcV(pLeaf1->m_translate, U_CA_C_TRANS);
	  pLeaf1->rotate(psis[i]);

	  pLeaf2 = new CLeaf(ROTAMER::UP, 0, D_PHI_AXIS, bbn_type, 2*i+2);
	  VcV(pLeaf2->m_translate, DU_C_CA_TRANS);
	  if (bbn_type == BBP)
	    pLeaf2->rotate(-M_PI/3.0);
	  else
	    pLeaf2->rotate(phis[i+1]);
	}

      m_links.push_back(pLeaf1);
      m_links.push_back(pLeaf2);

      last->setNext(pLeaf1);
      pLeaf1->setPrev(last);
      pLeaf1->setNext(pLeaf2);
      pLeaf2->setPrev(pLeaf1);
      last = pLeaf2;
    }
  last->setNext(NULL);

  m_frames.resize(getLength());
  computeTorsionEnergy();
}

REAL CChain::computeSolventRefEnergy()
{
  REAL sum = 0.0;

 for (int i = 0; i < getLength(); i++)
   {
     CLeaf * pLeaf = m_links[i];
     int type = pLeaf->getType();
     if (type == GLY)
	    continue;
     
     int size = SIDECHAIN::m_aalist[type]->m_size;
     for (int j = 0; j < size; j++)
       sum += deltaG_ref[SIDECHAIN::m_aalist[type]->m_aTypes[j]];
    }

 return sum;
}

void CChain::computeLinkCoords()
{
  vector<FRAME>::iterator it = m_frames.begin();
  Videntity(it->first.v);
  Midentity(it->second.m);

  REAL rot[3][3];
  REAL Rtemp[3][3];
  REAL trans[3];

  const CLeaf * pLeaf = getLink(0);
  
  McM(rot, pLeaf->m_rotate);
  VcV(trans, pLeaf->m_translate);

  while(pLeaf->getNext() != NULL)
    {
      it++;
      VcV(it->first.v, trans);
      McM(it->second.m, rot);

      pLeaf = pLeaf->getNext();

      MxVpV(trans, rot, pLeaf->m_translate, trans);
      
      McM(Rtemp, rot);
      MxM(rot, Rtemp,  pLeaf->m_rotate);
    }
}

void CChain::storeAngsStyle(const char * fname)
{
  ofstream fout(fname);

  if (!fout.is_open())
    {
      cout << " Could not open file " << fname << " for output!!!" << endl;
      return;
    }
  
  for (int j = 1; j < getLength(); j += 2)
  {
    fout << AA_NAMES[m_links[j]->getType()] << " ";
    REAL a = m_links[j-1]->getAngle() * 180/M_PI;
    while (a > 180.0)
      a -= 360.0;
    while (a < -180)
      a += 360;

    fout << a << " ";
    a = m_links[j]->getAngle() * 180/M_PI;
    while (a > 180.0)
      a -= 360.0;
    while (a < -180)
      a += 360;

    fout << a << " ";

    fout << m_links[j]->getRotIndex() << endl;
  }
}

// Store the chain in a PDB style file.
void CChain::storeCoordinates(const char * fname)
{
  computeLinkCoords();
  vector<AA> aas;
  AA aa;
  ATOM_ at;

  vector<FRAME>::const_iterator it = getFrames().begin();

  int i = 0, j = 0;
  for (j = 0; j < getLength(); j++)
    {
      if (m_links[j]->getType() == GLY)
	{
	  it++;
	  continue;
	}
    
      int size = SIDECHAIN::m_aalist[m_links[j]->getType()]->m_size;
      for (int k = 0; k < size; k++)
	{	
	  REAL vec[3];

	  if (j == 0) 
	    {
	      if (k == 0)
		aa.type = m_links[j+1]->getType();
	    }
	  else if (j % 2 == 0 && j != getLength() -1 && k == 2)
	      {
		aas.push_back(aa);
		aa.atoms.clear();
		aa.type = m_links[j+1]->getType();
	      }

	  MxVpV(vec, it->second.m, 
		m_links[j]->getPositions()[k], it->first.v);

	  strcpy(at.name, SIDECHAIN::m_aalist[m_links[j]->getType()]->m_aNames[k]);

	  at.pos[0] = vec[0];
	  at.pos[1] = vec[1];
	  at.pos[2] = vec[2];
	  aa.atoms.push_back(at);
	}
      
      it++;
    }
  aas.push_back(aa);
  
  writeToPDB(fname, aas);
}

// Compute the torsion energy, which is a function of the backbone angles
REAL CChain::computeTorsionEnergy()
{
  m_torsionE = 0.0;
  for (int i = 1; i < getLength() - 2; i++)
    m_torsionE += m_links[i]->getTorsionE();

  return m_torsionE;
}

// Update the torsion energy term after the gicen change is applied
// to the backbone angles.
REAL CChain::updateTorsionEnergy(const vector<ANGLE_CHANGE> & angles)
{
  m_undoTorsionE = m_torsionE;

  for (int i = 0; i < ((int)angles.size())-1; i++)
    m_torsionE += m_links[angles[i].m_index]->getTorsionChange();
    
  return m_torsionE;
}

//Compute the absolute positions of all Ca atoms in the protein chain.
void CChain::getCaPositions(POSITIONS & pos)
{
  int length = getLength() / 2;
  REAL Calphas[length][3];
  
  computeCaPositions(Calphas);

  if (pos.size() != length)
    pos.resize(length, vector<REAL>(3));

  for (int i = 0; i < length; i++)
    for (int j = 0; j < 3; j++)
      pos[i][j] = Calphas[i][j];
}

void CChain::computeCaPositions(REAL Calphas[][3])
{
  REAL rot[3][3];
  REAL Rtemp[3][3];
  REAL trans[3];
 
  const CLeaf * pLeaf = getLink(0);
  VcV(Calphas[0], pLeaf->getPositions()[4]);

  McM(rot, pLeaf->m_rotate);
  VcV(trans, pLeaf->m_translate);

  int i = 0;
  while(pLeaf->getNext()->getNext() != NULL)
    {
      pLeaf = pLeaf->getNext();

      if (pLeaf->getIndex() % 2 == 0)
	{
	  i++;
	  MxVpV(Calphas[i], rot, 
		pLeaf->getPositions()[4], trans);
	}

      MxVpV(trans, rot, pLeaf->m_translate, trans);      
      McM(Rtemp, rot);
      MxM(rot, Rtemp,  pLeaf->m_rotate);
    }
}

