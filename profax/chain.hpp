#ifndef _CHAIN_H
#define _CHAIN_H

#include <vector>
#include "pdb.hpp"
#include "leaf.hpp"

#define ANGS_FILE 1

struct REAL3 {
  REAL v[3];
};

struct REAL33 {
  REAL m[3][3];
};

typedef pair<REAL3, REAL33> FRAME;

class CSkipList;

class CChain
{
public:
  CChain(const char * fname, int type);

  void computeLinkCoords();
  void storeCoordinates(const char * fname);

  void storeAngsStyle(const char * fname);
  REAL computeTorsionEnergy();
  REAL updateTorsionEnergy(const vector<ANGLE_CHANGE> & angles);

  // Undo the change to the torsion energy caused by the latest move.
  void undoTorsionChange()
  { m_torsionE = m_undoTorsionE; }

  void getCaPositions(POSITIONS & pos);
  
  int getLength() const
  { return m_links.size(); }
  const vector<FRAME> & getFrames() const
  { return m_frames; }
  const vector<CLeaf*> & getLinks() const
  { return m_links; }
  CLeaf* getLink(int i) const
  { return m_links[i]; }
  REAL getTorsionEnergy() const
  { return m_torsionE; }
  REAL getSolventE() const
  { return m_solventE; }

private:
  void create(const vector<int> AAtypes, const vector<REAL> & phis, 
	      const vector<REAL> & psis, const vector<int> & rots);
  void load_angs(const char * fname, vector<int> & AAtypes, vector<REAL> & phis,
		 vector<REAL> & psis, vector<int> & rots);

  REAL computeSolventRefEnergy();
  void computeCaPositions(REAL Calphas[][3]);

  vector<FRAME> m_frames;
  vector<CLeaf*> m_links;

  REAL m_torsionE, m_undoTorsionE;
  REAL m_solventE;
};

#endif
