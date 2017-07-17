#ifndef _SLIST_H
#define _SLIST_H

#include "chain.hpp"
#include "node.hpp"
#include "leaf.hpp"
#include "pairtree.hpp"

class CSkiplist
{
public:
  // Create the hierarchy on top of the given chain.
  CSkiplist(CChain & chain) : 
    m_chain(chain), m_root(NULL), m_levels(0), m_time(0) 
    { createHierarchy(); }

  void clear();
  void makeMove(const vector<ANGLE_CHANGE> & angles, 
		const vector<ROTAMER_CHANGE> & rotamers);
  void undoLastMove();

  // Compute the absolute position of all atoms.
  // The first atom is at position (0,0,0).
  void computeAtomCenters() 
  { m_chain.computeLinkCoords(); }

  REAL computeDistance(int ind1, int i, int ind2, int j);
  REAL computeEnergy(const vector<ANGLE_CHANGE> & angles);

  bool findSelfClash();

  // Increment the timer.
  void incTime()
  { m_time++; }
 
  // Accessors
  const CNode * getRoot() const
  { return m_root; }
  const CChain & getChain() const
  { return m_chain; }
  int getLevels() const
  { return m_levels; }
  int getLength() const
  { return m_chain.getLength(); }
  REAL getEnergy() const
  { return m_energy; }
  CLeaf* getLink(int i)
  { return m_chain.getLink(i); }
  int getTime() const
  { return m_time; }

private:
  void createHierarchy();
  
  void clearUndo();
  void saveUndoNode(CNode * pNode)
  { m_undoNodes.push_back(pNode); }

  CNode * m_root;
  int m_levels;
  REAL m_energy;
  REAL m_undoEnergy;
  int m_time;

  vector<CNode*> m_undoNodes;

  CTerm * m_pairTree;
  CChain & m_chain;
};


#endif
