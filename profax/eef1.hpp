



#include <cmath>
#include <cassert>
#include <string.h>
#include "kernel.h"

#ifndef _EEF1_H
#define _EEF1_H
extern REAL CUTOFF_DISTANCE;
extern REAL CUTOFF_DISTANCE_2;
extern REAL CLASH_CUTOFF_DISTANCE;


extern char AA_NAMES[][5];
extern int NUM_ROTAMERS[];
extern int ROTAMER_START[];
extern float ROTAMER_VALUE[][4];

// for new g++ compilers use:
#include <ext/hash_map>
// for older versions use
//#include <hash_map>

#include "MatVec.h"



//#define T_0 298.15
//#define MAX_ROTAMER_SIZE 13 // size of ARG
//#define MAX_ROTAMER_PAIRS (MAX_ROTAMER_SIZE*(MAX_ROTAMER_SIZE-1))
////extern REAL CUTOFF_DISTANCE;
////extern REAL CUTOFF_DISTANCE_2;
////extern REAL CLASH_CUTOFF_DISTANCE;
//
//// The Boltzman constant
//#define Kb 0.00198576
//
////Atom types
//#define NUM_ATYPES 19
//#define NUM_HEAVY_TYPES 17
//
//#define  C    0
//#define  CR   1
//#define  CH1E 2
//#define  CH2E 3
//#define  CH3E 4
//#define  CR1E 5
//#define  NH1  6
//#define  NR   7
//#define  NH2  8
//#define  NH3  9
//#define  NC2  10
//#define  N    11
//#define  OH1  12
//#define  O    13
//#define  OC   14
//#define  S    15
//#define  SH1E 16
//#define  H    17
//#define  HC   18
//
//// Amino acid types
//#define NUMAA 25
//
//#define IND -1
//#define ARG 0
//#define ASN 1
//#define ASP 2
//#define CYS 3
//#define GLN 4
//#define GLU 5
//#define HIS 6
//#define ILE 7
//#define LEU 8
//#define LYS 9
//#define MET 10
//#define PHE 11
//#define PRO 12
//#define SER 13
//#define THR 14
//#define TRP 15
//#define TYR 16
//#define VAL 17
//#define GLY 18
//#define ALA 19
//#define NTR 20 // C-terminal cap
//#define CTR 21 // N-terminal cap
//#define BBN 22 // C - N backbone piece
//#define BBP 23 // C - N backbone piece for Proline
//#define BBG 24 // C - N backbone piece for Glycine
//
//extern char AA_NAMES[][5];
//extern int NUM_ROTAMERS[];
//extern int ROTAMER_START[];
//extern float ROTAMER_VALUE[][4];
//
extern REAL epsilon[NUM_ATYPES];
extern REAL SIGMA[NUM_ATYPES];
//
//// The factor used to determine when we have a steric clash
//#define SIGMA_FACTOR 0.5
//
extern const REAL DIELECTRIC;

int getAA(const char * aa);



inline
REAL compute_ES(REAL charge1, REAL charge2, REAL dist)
{
  // We assume the "dist" is the squared distance.
  return (DIELECTRIC * (charge1 * charge2)/(dist));
}

inline
bool isStericClash(int type1, int type2, REAL dist, bool is14)
{
  REAL sig1 = SIGMA[type1];
  REAL sig2 = SIGMA[type2];

  if (is14)
    {
      if (type1 < NH1)
	  sig1 = 1.9;

      if (type2 < NH1)
	  sig2 = 1.9;
    }

  REAL sig = SIGMA_FACTOR * (sig1 + sig2);

  return (dist < sig*sig);
}

inline
REAL compute_vdW(int type1, int type2, REAL dist, bool is14)
{
  // We assume "dist" is the distance squared.
  if (dist > CUTOFF_DISTANCE_2)
    return 0.0;

  REAL sig1 = SIGMA[type1];
  REAL eps1 = epsilon[type1];
  REAL sig2 = SIGMA[type2];
  REAL eps2 = epsilon[type2];

  if (is14)
    {
      if (type1 < NH1)
	{
	  sig1 = 1.9;
	  eps1 = -0.1;
	}

      if (type2 < NH1)
	{
	  sig2 = 1.9;
	  eps2 = -0.1;
	}
    }

  REAL eps = sqrt(eps1*eps2);
  REAL sig = (sig1 + sig2);
 
  REAL rat2 = (sig*sig)/dist;
  REAL rat6 = rat2*rat2*rat2;
		   
  REAL res = (eps*(rat6*rat6 - 2*rat6));   
  return res;
}


// Atom volumes
extern REAL volume[NUM_ATYPES];
extern REAL deltaG_ref[NUM_ATYPES];
extern REAL deltaG_free[NUM_HEAVY_TYPES];
extern REAL deltaH_ref[NUM_HEAVY_TYPES];
extern REAL deltaCp_ref[NUM_HEAVY_TYPES];

inline
REAL compute_deltaG_ref(REAL T, int aType)
{
  assert(aType < NUM_HEAVY_TYPES && aType >= 0);

  REAL deltaS_ref = (deltaH_ref[aType] - deltaG_ref[aType]) / T_0;
  
  REAL res = deltaG_ref[aType] - deltaS_ref*(T - T_0) - 
    deltaCp_ref[aType]*T*log(T/T_0) + deltaCp_ref[aType]*(T-T_0);

  return res;
}

inline
REAL getLambda(int AAtype, int index, int aType)
{
  if (aType == NH3 || aType == NC2 || aType == OC)
    return 6.0;

  if ((AAtype == ARG && index >= 2) ||
      (AAtype == LYS && index >= 3) ||
      (AAtype == ASP) ||
      (AAtype == GLU && index >= 1) ||
      (AAtype == NTR) ||
      (AAtype = CTR))
    return 6.00;
  else 
    return 3.50;
}

extern const REAL SOLVATION_K;
inline
REAL computeSolventEffect(int aType1, int aType2, REAL dist, REAL lambda1,
			  REAL lambda2, bool is14)
{
  if (dist > CUTOFF_DISTANCE_2)
    return 0.0;

  REAL t1 = 0.0, t2 = 0.0;

  if (aType1 < NUM_HEAVY_TYPES && aType2 < NUM_HEAVY_TYPES)
  {
    REAL sig1 = SIGMA[aType1], sig2 = SIGMA[aType2];

  if (is14)
    {
      if (aType1 < NH1)
	sig1 = 1.9;

      if (aType2 < NH1)
	sig2 = 1.9;
    }

    // We assume "dist" is the squared distance.
    REAL d = sqrt(dist);

    REAL X12 = (d - sig1)/lambda1;
    t1 = SOLVATION_K * deltaG_free[aType1] * exp(-X12*X12) * volume[aType2] 
      / (lambda1 * dist);

    REAL X21 = (d - sig2)/lambda2;
    t2 = SOLVATION_K * deltaG_free[aType2] * exp(-X21*X21) * volume[aType1] 
      / (lambda2 * dist);
  }

  return -(t2 + t1);
}

// dihedral constants
#define NUM_DIHEDRALS 6

extern const REAL E_0[NUM_DIHEDRALS];
extern const REAL theta_0[NUM_DIHEDRALS];
extern const REAL dihed_n[NUM_DIHEDRALS];

inline
REAL compute_dihedral(REAL angle, int type)
{
  assert(type >= 0 && type < NUM_DIHEDRALS);

  return E_0[type] * (1 + cos(dihed_n[type]*angle - theta_0[type]));
}

class CBV;

// A structure to hold rotamer information that can be switched quickly
// into a leaf.
struct ROTAMER
{
	ROTAMER() : m_energy(0.0), m_bv(NULL)
	    {
	        m_positions = new REAL[0][3];
	    }

    /*ROTAMER ():ROTAMER(0) {
        
    }*/

    ROTAMER(int size) : m_energy(0.0), m_bv(NULL)
    {
        m_positions = new REAL[size][3];
    }
    
    enum ROT_TYPE {UP = 1, DOWN = -1};
    REAL (*m_positions)[3];
    CBV * m_bv;
    REAL m_energy;
};

REAL compute_rotamer_energy(int type, const ROTAMER & rot, int index);

// A structur to hold all sidechain information
struct SIDECHAIN
{
  SIDECHAIN(int size, int nRotamers, int nGroups, int nChis) : 
    m_size(size), m_nRotamers(nRotamers), m_nGroups(nGroups), 
    m_nChis(nChis)
  {
    m_Urotamers = new ROTAMER[nRotamers];
    //m_Urotamers = new ROTAMER(size);
    m_Drotamers = new ROTAMER [nRotamers];
    //m_Drotamers = new ROTAMER(size);

      
      for(int i = 0; i< nRotamers; i++) {
          m_Drotamers[i] = ROTAMER(size);
          m_Urotamers[i] = ROTAMER(size);
      }
      
    m_groups = new int[nGroups];
    if (nChis > 0)
      m_chiTypes = new int[nChis];
    else 
      m_chiTypes = NULL;
  }

  const ROTAMER & getRotamer(int index, ROTAMER::ROT_TYPE rType)
  {
    if (rType == ROTAMER::UP)
      return m_Urotamers[index];
    else if  (rType == ROTAMER::DOWN)
      return m_Drotamers[index];
    else
      assert(false);
  } 

  static void createSidechains(const char * dir);
  static SIDECHAIN * m_aalist[NUMAA];

  ROTAMER * m_Urotamers;
  ROTAMER * m_Drotamers;
  int m_nRotamers;
  int m_size;
  int m_nGroups;
  int * m_groups;
  REAL * m_charges;
  int * m_aTypes;
  int m_nChis;
  int * m_chiTypes;
  char ** m_aNames;
};

// Axis of rotatin for PHI and PSI angles
extern const REAL U_PHI_AXIS[3];
extern const REAL D_PHI_AXIS[3];

extern const REAL U_PSI_AXIS[3];
extern const REAL D_PSI_AXIS[3];

// Translation between origins of the different kinds of links.
extern const REAL N_C_NTR_TRANS[3];

extern const REAL U_CA_C_TRANS[3];
extern const REAL D_CA_C_TRANS[3];

extern const REAL UD_C_CA_TRANS[3];
extern const REAL DU_C_CA_TRANS[3];

// Exclusion list
inline
int to_index(int diff, int t1, int i, int t2, int j)
{
  return (diff*1024*1024 + t1*1024*32 + i*1024 + t2*32 + j);
}

// for new g++ versions use:
typedef __gnu_cxx::hash_map<int, int> EXCLUSIONS;
// for older versions use
//typedef hash_map<int, int> EXCLUSIONS;

extern EXCLUSIONS exclusion_list;
#define EXCLUDED 1
#define PAIR1_4  2
#define NOT_EXCLUDED 0

inline
int isExcluded(int diff, int type1, int id1, int type2, int id2)
{
  if (diff > 3) 
    return NOT_EXCLUDED;
  else if (diff == 0 && id1 >= id2)
    return EXCLUDED;

  EXCLUSIONS::const_iterator ex = exclusion_list.find(to_index(diff,type1,id1,type2,id2));

  if (ex == exclusion_list.end())
    return NOT_EXCLUDED;
  else if (ex->second < 3)
    return EXCLUDED;
  else if (ex->second == 3)
    return PAIR1_4;
  else
    assert(false);
}

void Initialize(const char * dir);

#endif
