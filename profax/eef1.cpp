#include <fstream>
#include <utility>
#include <iostream>

#include "eef1.hpp"
#include "rss.hpp"
#include "spheres.hpp"
#include "pairtree.hpp"

using namespace std;


#ifndef _EEF1_C_
#define _EEF1_C_

#define ROTAMER_FILE "rotamer_coords"

char AA_NAMES[NUMAA][5]={"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
			 "HIS", "ILE", "LEU", "LYS", "MET", "PHE", 
			 "PRO", "SER", "THR", "TRP", "TYR", "VAL", 
			 "GLY", "ALA", "NTR", "CTR", "BBN", "BBP",
			 "BBG"};
int NUM_ROTAMERS[NUMAA] = {81, 9, 9, 3, 27, 27, 6, 9, 9, 81, 27, 
			   6, 2, 3, 3, 9, 6, 3, 1, 1, 1, 1, 1, 
			   1, 1};
int ROTAMER_START[NUMAA] = {0, 81, 90, 99, 102, 129, 156, 162, 171, 180, 
			    261, 288, 294, 296, 299, 302, 311, 317, 320,
			    320, 320, 320, 320, 320, 320};


float ROTAMER_VALUE[][4]=
#include "rotamers.h"

// van der Waals constants
REAL epsilon[NUM_ATYPES] = {-0.1200, -0.1200, -0.0486, -0.1142, -0.1811, 
			    -0.1200, -0.2384, -0.2384, -0.2384, -0.2384,
			    -0.2384, -0.2384, -0.1591, -0.1591, -0.6469,
			    -0.0430, -0.0430, -0.0498, -0.0498};

REAL SIGMA[NUM_ATYPES] = {2.100, 2.100, 2.365, 2.235, 2.165, 2.100, 
			  1.6000, 1.6000, 1.6000, 1.6000, 1.6000, 1.6000,
			  1.6000, 1.6000, 1.6000,  1.890, 1.890, 0.8000, 
			  0.6000};

REAL CUTOFF_DISTANCE = 9.0;
REAL CUTOFF_DISTANCE_2 = CUTOFF_DISTANCE*CUTOFF_DISTANCE;
REAL CLASH_CUTOFF_DISTANCE = SIGMA[2];

// Atom volumes
REAL volume[NUM_ATYPES] = {14.7, 8.3, 23.7, 22.4, 30.0, 18.4, 4.4, 4.4, 
			   11.2, 11.2, 11.2, 0.0, 10.8, 10.8, 10.8, 14.7, 21.4,
			   0.0, 0.0};

// Solvation energy parameters

// Atom \Delta G^{ref}
REAL deltaG_ref[NUM_ATYPES] = {0.000, -0.890, -0.187, 0.372, 1.089,
			       0.057, -5.950, -3.820, -5.450, -20.000,
			       -10.000, -1.000, -5.920, -5.330, 
			       -10.000, -3.240, -2.050, 0.0, 0.0};

// Atom \Delta G^{free}
REAL deltaG_free[NUM_HEAVY_TYPES] = {0.00, -1.40, -0.25, 0.52, 1.50, 0.08,
				-8.90, -4.00, -7.80, -20.00, -10.00,
				-1.55, -6.70, -5.85, -10.00, -4.10, 
				-2.70}; 

// Atom \Delta H^{ref}
REAL deltaH_ref[NUM_HEAVY_TYPES] = {0.000, 2.220, 0.876, -0.610, -1.779,
			       -0.973, -9.059, -4.654, -9.028, -25.000,
			       -12.000, -1.250, -9.264, -5.787, -12.000,
			       -4.475, -4.475};

// Atom \Delta Cp^{ref}
REAL deltaCp_ref[NUM_HEAVY_TYPES] = {0.00, 6.90, 0.00, 18.60, 35.60, 6.90,
				-8.80, -8.80, -7.00, -18.00, -7.00, 8.80,
				-11.20, -8.80, -9.40, -39.90, -39.90};

const REAL SOLVATION_K = 2.0 / (4.0 * M_PI * sqrt(M_PI));

// dihedral constants
const REAL E_0[NUM_DIHEDRALS] = {0.3, 0, 1.6, 0.3, 0.0, 1.2}; 
const REAL theta_0[NUM_DIHEDRALS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const REAL dihed_n[NUM_DIHEDRALS] = {3, 3, 3, 3, 3, 2};

#endif

// Convert AA name to a numerical code.
int getAA(const char * aa)
{
  if (strlen(aa) == 3)
    {
      if (strcmp(aa, "ALA") == 0)
	return (ALA);
      else if (strcmp(aa, "ARG") == 0)
	return (ARG);
      else if (strcmp(aa, "ASN") == 0)
	return (ASN);
      else if (strcmp(aa, "ASP") == 0)
	return (ASP);
      else if (strcmp(aa, "CYS") == 0)
	return (CYS);
      else if (strcmp(aa, "GLN") == 0)
	return (GLN);
      else if (strcmp(aa, "GLU") == 0)
	return (GLU);
      else if (strcmp(aa, "GLY") == 0)
	return (GLY);
      else if (strcmp(aa, "HIS") == 0)
	return (HIS);
      else if (strcmp(aa, "ILE") == 0)
	return (ILE);
      else if (strcmp(aa, "LEU") == 0)
	return (LEU);
      else if (strcmp(aa, "LYS") == 0)
	return (LYS);
      else if (strcmp(aa, "MET") == 0)
	return (MET);
      else if (strcmp(aa, "PHE") == 0)
	return (PHE);
      else if (strcmp(aa, "PRO") == 0)
	return (PRO);
      else if (strcmp(aa, "SER") == 0)
	return (SER);
      else if (strcmp(aa, "THR") == 0)
	return (THR);
      else if (strcmp(aa, "TRP") == 0)
	return (TRP);
      else if (strcmp(aa, "TYR") == 0)
	return (TYR);
      else if (strcmp(aa, "VAL") == 0)
	return (VAL);
      else if (strcmp(aa, "CTR") == 0)
	return (CTR);
      else if (strcmp(aa, "NTR") == 0)
	return (NTR);
      else if (strcmp(aa, "BBN") == 0)
	return (BBN);
      else if (strcmp(aa, "BBG") == 0)
	return (BBG);
      else if (strcmp(aa, "BBP") == 0)
	return (BBP);
      else
	{
	  cout << "bad AA name " << aa << endl;
	  exit(0);
	}
    }
  else if (strlen(aa) == 1)
    {
      if (aa[0] == 'A')
	return ALA;
      else if (aa[0] == 'C')
	return CYS;
      else if (aa[0] == 'D')
	return ASP;
      else if (aa[0] == 'E')
	return GLU;
      else if (aa[0] == 'F')
	return PHE;
      else if (aa[0] == 'G')
	return GLY;
      else if (aa[0] == 'H')
	return HIS;
      else if (aa[0] == 'I')
	return ILE;
      else if (aa[0] == 'K')
	return LYS;
      else if (aa[0] == 'L')
	return LEU;
      else if (aa[0] == 'M')
	return MET;
      else if (aa[0] == 'N')
	return ASN;
      else if (aa[0] == 'P')
	return PRO;
      else if (aa[0] == 'Q')
	return GLN;
      else if (aa[0] == 'R')
	return ARG;
      else if (aa[0] == 'S')
	return SER;
      else if (aa[0] == 'T')
	return THR;
      else if (aa[0] == 'V')
	return VAL;
      else if (aa[0] == 'W')
	return TRP;
      else if (aa[0] == 'Y')
	return TYR;
      else
	{
	  cout << "bad AA name " << aa << endl;
	  exit(0);
	}
    }
  else
    {
      cout << "bad AA name " << aa << endl;
      exit(0);
    }
}

// Vacuum Dielectric constant
const REAL DIELECTRIC = 332.05382;


// Axes of rotation for PHI and PSI angles
const REAL D_PHI_AXIS[3] = {0.0158, 0.0000, 0.9999};
const REAL U_PHI_AXIS[3] = {0.9766, 0.0000, 0.2150};

const REAL D_PSI_AXIS[3] = {0.9353, 0.0000, 0.3537};
const REAL U_PSI_AXIS[3] = {0.1597, 0.0000, 0.9872};

// Translation between origins of the different kinds of links.
const REAL N_C_NTR_TRANS[3] = {0.0000, 0.0000, 1.5000};

const REAL D_CA_C_TRANS[3] = {1.4120, 0.0000, 0.5340};
const REAL U_CA_C_TRANS[3] = {0.2410, 0.0000, 1.4900};

const REAL UD_C_CA_TRANS[3] = {2.3470, 0.0000, -0.6600};
const REAL DU_C_CA_TRANS[3] = {-1.1150, 0.0000, 2.1680};
		  
SIDECHAIN * SIDECHAIN::m_aalist[NUMAA];
#include "AA.hpp"

// Create an array of all Amino acids each with all its rotamers.
// Each rotamer is saved in two configurations: UP and DOWN to mimic 
// the effect of the 180 deg Omega angle.
void SIDECHAIN::createSidechains(const char * dir)
{
  char name[100];
  sprintf(name, "%s/%s", dir ,ROTAMER_FILE);
  ifstream fin(name);
     
  if (!fin.is_open())
    {
      cout << "Could not open file of rotamers!!!" << endl;
      exit(0);
    }

  // Alanine
  SIDECHAIN * pSC = new SIDECHAIN(ALA_size, NUM_ROTAMERS[ALA], ALA_nGroups,
				  ALA_nChis); 
  pSC->m_charges = ALA_charges;
  pSC->m_groups = ALA_groups;
  pSC->m_aTypes = ALA_aTypes;
  pSC->m_chiTypes = ALA_chiTypes;
  pSC->m_aNames = ALA_aNames;
  m_aalist[ALA] = pSC;

  // Arginine
  pSC = new SIDECHAIN(ARG_size, NUM_ROTAMERS[ARG], ARG_nGroups, ARG_nChis); 
  pSC->m_charges = ARG_charges;
  pSC->m_groups = ARG_groups;
  pSC->m_aTypes = ARG_aTypes;
  pSC->m_chiTypes = ARG_chiTypes;
  pSC->m_aNames = ARG_aNames;
  m_aalist[ARG] = pSC;

  // Aspartic acid
  pSC = new SIDECHAIN(ASN_size, NUM_ROTAMERS[ASN], ASN_nGroups, ASN_nChis); 
  pSC->m_charges = ASN_charges;
  pSC->m_groups = ASN_groups;
  pSC->m_aTypes = ASN_aTypes;
  pSC->m_chiTypes = ASN_chiTypes;
  pSC->m_aNames = ASN_aNames;
  m_aalist[ASN] = pSC;

  // Aspartine
  pSC = new SIDECHAIN(ASP_size, NUM_ROTAMERS[ASP], ASP_nGroups, ASP_nChis); 
  pSC->m_charges = ASP_charges;
  pSC->m_groups = ASP_groups;
  pSC->m_aTypes = ASP_aTypes;
  pSC->m_chiTypes = ASP_chiTypes;
  pSC->m_aNames = ASP_aNames;
  m_aalist[ASP] = pSC;

  // Cystine
  pSC = new SIDECHAIN(CYS_size, NUM_ROTAMERS[CYS], CYS_nGroups, CYS_nChis); 
  pSC->m_charges = CYS_charges;
  pSC->m_groups = CYS_groups;
  pSC->m_aTypes = CYS_aTypes;
  pSC->m_chiTypes = CYS_chiTypes;
  pSC->m_aNames = CYS_aNames;
  m_aalist[CYS] = pSC;

  // Glutamic acid
  pSC = new SIDECHAIN(GLN_size, NUM_ROTAMERS[GLN], GLN_nGroups, GLN_nChis); 
  pSC->m_charges = GLN_charges;
  pSC->m_groups = GLN_groups;
  pSC->m_aTypes = GLN_aTypes;
  pSC->m_chiTypes = GLN_chiTypes;
  pSC->m_aNames = GLN_aNames;
  m_aalist[GLN] = pSC;

  // Glutamine
  pSC = new SIDECHAIN(GLU_size, NUM_ROTAMERS[GLU], GLU_nGroups, GLU_nChis); 
  pSC->m_charges = GLU_charges;
  pSC->m_groups = GLU_groups;
  pSC->m_aTypes = GLU_aTypes;
  pSC->m_chiTypes = GLU_chiTypes;
  pSC->m_aNames = GLU_aNames;
  m_aalist[GLU] = pSC;

  // Glysine
  pSC = new SIDECHAIN(GLY_size, NUM_ROTAMERS[GLY], GLY_nGroups, GLY_nChis); 
  pSC->m_charges = GLY_charges;
  pSC->m_groups = GLY_groups;
  pSC->m_aTypes = GLY_aTypes;
  pSC->m_chiTypes = GLY_chiTypes;
  pSC->m_aNames = GLY_aNames;
  m_aalist[GLY] = pSC;

  // Histidine
  pSC = new SIDECHAIN(HIS_size, NUM_ROTAMERS[HIS], HIS_nGroups, HIS_nChis); 
  pSC->m_charges = HIS_charges;
  pSC->m_groups = HIS_groups;
  pSC->m_aTypes = HIS_aTypes;
  pSC->m_chiTypes = HIS_chiTypes;
  pSC->m_aNames = HIS_aNames;
  m_aalist[HIS] = pSC;

  // Isoleucine
  pSC = new SIDECHAIN(ILE_size, NUM_ROTAMERS[ILE], ILE_nGroups, ILE_nChis); 
  pSC->m_charges = ILE_charges;
  pSC->m_groups = ILE_groups;
  pSC->m_aTypes = ILE_aTypes;
  pSC->m_chiTypes = ILE_chiTypes;
  pSC->m_aNames = ILE_aNames;
  m_aalist[ILE] = pSC;

  // Leucine
  pSC = new SIDECHAIN(LEU_size, NUM_ROTAMERS[LEU], LEU_nGroups, LEU_nChis); 
  pSC->m_charges = LEU_charges;
  pSC->m_groups = LEU_groups;
  pSC->m_aTypes = LEU_aTypes;
  pSC->m_chiTypes = LEU_chiTypes;
  pSC->m_aNames = LEU_aNames;
  m_aalist[LEU] = pSC;

  // Lysine
  pSC = new SIDECHAIN(LYS_size, NUM_ROTAMERS[LYS], LYS_nGroups, LYS_nChis); 
  pSC->m_charges = LYS_charges;
  pSC->m_groups = LYS_groups;
  pSC->m_aTypes = LYS_aTypes;
  pSC->m_chiTypes = LYS_chiTypes;
  pSC->m_aNames = LYS_aNames;
  m_aalist[LYS] = pSC;

  // Methionine
  pSC = new SIDECHAIN(MET_size, NUM_ROTAMERS[MET], MET_nGroups, MET_nChis); 
  pSC->m_charges = MET_charges;
  pSC->m_groups = MET_groups;
  pSC->m_aTypes = MET_aTypes;
  pSC->m_chiTypes = MET_chiTypes;
  pSC->m_aNames = MET_aNames;
  m_aalist[MET] = pSC;

  // Phenylalanine
  pSC = new SIDECHAIN(PHE_size, NUM_ROTAMERS[PHE], PHE_nGroups, PHE_nChis); 
  pSC->m_charges = PHE_charges;
  pSC->m_groups = PHE_groups;
  pSC->m_aTypes = PHE_aTypes;
  pSC->m_chiTypes = PHE_chiTypes;
  pSC->m_aNames = PHE_aNames;
  m_aalist[PHE] = pSC;
 
  // Proline
  pSC = new SIDECHAIN(PRO_size, NUM_ROTAMERS[PRO], PRO_nGroups, PRO_nChis); 
  pSC->m_charges = PRO_charges;
  pSC->m_groups = PRO_groups;
  pSC->m_aTypes = PRO_aTypes;
  pSC->m_chiTypes = PRO_chiTypes;
  pSC->m_aNames = PRO_aNames;
  m_aalist[PRO] = pSC;
 
  // Serine
  pSC = new SIDECHAIN(SER_size, NUM_ROTAMERS[SER], SER_nGroups, SER_nChis); 
  pSC->m_charges = SER_charges;
  pSC->m_groups = SER_groups;
  pSC->m_aTypes = SER_aTypes;
  pSC->m_chiTypes = SER_chiTypes;
  pSC->m_aNames = SER_aNames;
  m_aalist[SER] = pSC;

  // Thrionine
  pSC = new SIDECHAIN(THR_size, NUM_ROTAMERS[THR], THR_nGroups, THR_nChis); 
  pSC->m_charges = THR_charges;
  pSC->m_groups = THR_groups;
  pSC->m_aTypes = THR_aTypes;
  pSC->m_chiTypes = THR_chiTypes;
  pSC->m_aNames = THR_aNames;
  m_aalist[THR] = pSC;
 
  // Triptophan
  pSC = new SIDECHAIN(TRP_size, NUM_ROTAMERS[TRP], TRP_nGroups, TRP_nChis); 
  pSC->m_charges = TRP_charges;
  pSC->m_groups = TRP_groups;
  pSC->m_aTypes = TRP_aTypes;
  pSC->m_chiTypes = TRP_chiTypes;
  pSC->m_aNames = TRP_aNames;
  m_aalist[TRP] = pSC;
 
  // Tyrosine
  pSC = new SIDECHAIN(TYR_size, NUM_ROTAMERS[TYR], TYR_nGroups, TYR_nChis); 
  pSC->m_charges = TYR_charges;
  pSC->m_groups = TYR_groups;
  pSC->m_aTypes = TYR_aTypes;
  pSC->m_chiTypes = TYR_chiTypes;
  pSC->m_aNames = TYR_aNames;
  m_aalist[TYR] = pSC;
 
  // Valine
  pSC = new SIDECHAIN(VAL_size, NUM_ROTAMERS[VAL], VAL_nGroups, VAL_nChis); 
  pSC->m_charges = VAL_charges;
  pSC->m_groups = VAL_groups;
  pSC->m_aTypes = VAL_aTypes;
  pSC->m_chiTypes = VAL_chiTypes;
  pSC->m_aNames = VAL_aNames;
  m_aalist[VAL] = pSC;

  // C-terminal
  pSC = new SIDECHAIN(CTR_size, NUM_ROTAMERS[CTR], CTR_nGroups, 0); 
  pSC->m_charges = CTR_charges;
  pSC->m_groups = CTR_groups;
  pSC->m_aTypes = CTR_aTypes;
  pSC->m_aNames = CTR_aNames;
  m_aalist[CTR] = pSC;

  // N-terminal 
  pSC = new SIDECHAIN(NTR_size, NUM_ROTAMERS[NTR], NTR_nGroups, 0); 
  pSC->m_charges = NTR_charges;
  pSC->m_groups = NTR_groups;
  pSC->m_aTypes = NTR_aTypes;
  pSC->m_aNames = NTR_aNames;
  m_aalist[NTR] = pSC;

  // C - N backbone
  pSC = new SIDECHAIN(BBN_size, NUM_ROTAMERS[BBN], BBN_nGroups, 0); 
  pSC->m_charges = BBN_charges;
  pSC->m_groups = BBN_groups;
  pSC->m_aTypes = BBN_aTypes;
  pSC->m_aNames = BBN_aNames;
  m_aalist[BBN] = pSC;

  // C - N backbone for Proline.
  pSC = new SIDECHAIN(BBP_size, NUM_ROTAMERS[BBP], BBP_nGroups, 0); 
  pSC->m_charges = BBP_charges;
  pSC->m_groups = BBP_groups;
  pSC->m_aTypes = BBP_aTypes;
  pSC->m_aNames = BBP_aNames;
  m_aalist[BBP] = pSC;

  // C - N backbone for Glycine.
  pSC = new SIDECHAIN(BBG_size, NUM_ROTAMERS[BBG], BBG_nGroups, 0); 
  pSC->m_charges = BBG_charges;
  pSC->m_groups = BBG_groups;
  pSC->m_aTypes = BBG_aTypes;
  pSC->m_aNames = BBG_aNames;
  m_aalist[BBG] = pSC;

  // Load the coordinates of all rotamers of all AAs.
  char buf[100];
  for (int aa = 0; aa < NUMAA; aa++)
    {
      pSC = m_aalist[aa];
      for (int i = 0; i < NUM_ROTAMERS[aa]; i++)
	{
	  fin >> buf; // read separator line.

	  // Load the UP rotamer
	  for (int j = 0; j < pSC->m_size; j++)
	    {
	      fin >> pSC->m_Urotamers[i].m_positions[j][0];
	      fin >> pSC->m_Urotamers[i].m_positions[j][1];
	      fin >> pSC->m_Urotamers[i].m_positions[j][2];
	    }

	  // Compute the BV for this rotamer
#ifdef USE_RSS
	  pSC->m_Urotamers[i].m_bv = 
	    (CBV*) new CRss(pSC->m_Urotamers[i].m_positions,  pSC->m_size);
#else
	  pSC->m_Urotamers[i].m_bv = 
	    (CBV*) new CSphere(pSC->m_Urotamers[i].m_positions,  pSC->m_size);
#endif

	  // compute the energy of the rotamer
	  pSC->m_Urotamers[i].m_energy = 
	    compute_rotamer_energy(aa, pSC->m_Urotamers[i], i);

	  fin >> buf; // read separator line.

	  // Load the DOWN rotamer
	  for (int j = 0; j < pSC->m_size; j++)
	    {
	      fin >> pSC->m_Drotamers[i].m_positions[j][0];
	      fin >> pSC->m_Drotamers[i].m_positions[j][1];
	      fin >> pSC->m_Drotamers[i].m_positions[j][2];
	    }

	  // Compute the RSS for this rotamer.
#ifdef USE_RSS
	  pSC->m_Drotamers[i].m_bv = 
	    (CBV*) new CRss(pSC->m_Drotamers[i].m_positions,  pSC->m_size);	 
#else
	  pSC->m_Drotamers[i].m_bv = 
	    (CBV*) new CSphere(pSC->m_Drotamers[i].m_positions,  pSC->m_size);	 
#endif

	  // The energy of the D-rotamer is the same as the U-rotamer
	  pSC->m_Drotamers[i].m_energy = pSC->m_Urotamers[i].m_energy;
	}
    }
}

REAL compute_rotamer_energy(int type, const ROTAMER & rot, int index)
{
  // Glycine has no sidechain.
  if (type == GLY)
    return 0.0;

  REAL dists[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];

  int size = SIDECHAIN::m_aalist[type]->m_size;

  // Compute the full distances matrix
  REAL cen[3], dist[3];
  for(int i = 0; i < size; i++)
    for (int j = i+1; j < size; j++)
      {
	VmV(dist, rot.m_positions[i], rot.m_positions[j]);
	dists[i][j] = Vlength2(dist);
	assert(dists[i][j] > 0.0);
      }

  REAL sum = 0.0;

  // Compute all vdW terms.
  sum += CTerm::computeVdW(type, type, dists, 0);

  // Compute all elctrostatic terms
  sum += CTerm::computeElectrostatics(type, type, dists, 0);

  // Compute all solvation terms.
  sum += CTerm::computeSolvation(type, type, dists, 0);

  // Compute dihedral contribution
  int nChis = SIDECHAIN::m_aalist[type]->m_nChis;
  REAL s = 0.0;
  for (int a = 0; a < nChis; a++)
    {
      REAL angle = ROTAMER_VALUE[ROTAMER_START[type] + index][a];
      int dType = SIDECHAIN::m_aalist[type]->m_chiTypes[a];

      s += compute_dihedral(angle * M_PI/180, dType);
    }

  sum += s;

  return sum;
}

EXCLUSIONS exclusion_list;

// Create the list of excluded interactions (all atoms 2 bonds apart or less)
// 1-4 atoms (3 bonds apart) are in the list but are not excluded. In some cases
// they are treated differently.
void createExclusionList(const char * dir)
{
  char buf[250];
  sprintf(buf, "%s/exclusions", dir);
  ifstream fin(buf);
  if (!fin.is_open())
    {
      cout << "Could not open exclusions file: " << buf << endl;
      exit(0);
    }

  int diff, id1, id2, bonds;
  char name1[4], name2[4];
  while(!fin.eof())
    {
      fin.getline(buf, 250);

      if (buf[0] != '#')
	{
	  int res = sscanf(buf, "%d %s %d %s %d %d", &diff, name1, &id1, 
			   name2, &id2, &bonds);
	  assert(res == 6); 
	  exclusion_list.insert(make_pair(to_index(diff, getAA(name1), id1, 
						   getAA(name2),
						   id2), bonds));
	}
    }
}

// Initialize the exclusion list and all the sidechain data structures.
void Initialize(const char * dir)
{
  createExclusionList(dir);
  SIDECHAIN::createSidechains(dir); 
}

