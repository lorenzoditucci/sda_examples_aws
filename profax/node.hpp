#ifndef _NODE_H
#define _NODE_H

#include "eef1.hpp"
#include "bv.hpp"

#include <vector>

using namespace std;

class CSkiplist;
class CTerm;

class CNode
{
public:
  // Construct an empty node withthe given index.
  CNode(int index) :  m_next(NULL), m_prev(NULL), m_up(NULL), 
		      m_down(NULL), m_index(index), m_stamp(0), 
		      m_level(0), m_nc(BOTH), m_pSL(NULL), m_bv(NULL)
  { Midentity(m_rotate); Videntity(m_translate); } 
  
  CNode(CNode * down, CSkiplist * pSL); 

  // The types of changes that a node can undergo.
  enum NODE_CHANGE {NONE = 0, TRANSFORM = 1, BOX = 2, BOTH = 3};
  
  void updateBV();

  // Undo the changes to the BV casued by the last move.
  void undoBV()
  { m_bv->undoBV(); }

  void transformToNext();
  virtual void computePairEnergy(CNode * pNode, const REAL rot[3][3],
				 const REAL trans[3], CTerm * term,
				 bool bSeparated);
  virtual void computeSelfEnergy(CTerm * term);
  virtual bool findPairClash(CNode * pNode, const REAL rot[3][3],
			     const REAL trans[3], bool bSeparated);
  virtual bool findSelfClash();

  virtual void undo();

  // Modifiers
  void setBV(CBV * bv)
  { m_bv = bv; }
  void setNext(CNode * next)
  { m_next = next; }
  void setPrev(CNode * prev)
  { m_prev = prev; }
  void setUp(CNode * up)
  { m_up = up; }
  void setDown(CNode * down)
  { m_down = down; }
  void setLevel(int l)
  { m_level = l; }
  void setStamp(int stamp, NODE_CHANGE nc)
  { m_stamp = stamp; m_nc = nc; }

  // Accessors
  CBV * getBV() const
  { return m_bv; }
  CNode * getNext() const
  { return m_next; }
  CNode * getPrev() const
  { return m_prev; }
  CNode * goUp() const
  { return m_up; }
  CNode * goDown() const
  { return m_down; }
  CNode * getParent() const
  { 
    if (goUp())
	return goUp();
    else if (getPrev())
      return getPrev()->goUp();
    else
      return NULL;
  }
  int getLevel() const
  { return m_level; }
  int getIndex() const
  { return m_index; }
  int getStamp() const
  { return m_stamp; }
  virtual bool isLeaf() const
  { return false; }
  bool isSeparator(int time)
  { return((m_stamp >= time) && (m_nc & TRANSFORM)); }
  bool isAffected(int time)
  { return ((m_stamp >= time) && (m_nc & BOX)); }

  const REAL getBVVolume() const
  { assert(m_bv); return m_bv->getVolume(); }
 
  friend ostream & operator<<(ostream & out, CNode & b);

  static vector<CTerm*> m_undoPairs;

  REAL m_rotate[3][3];
  REAL m_translate[3];
  CBV * m_bv;
  CSkiplist * m_pSL;

protected:
  CNode * m_next;
  NODE_CHANGE m_nc;

  // undo variables
  REAL m_undoRotate[3][3];
  REAL m_undoTranslate[3];

private:
  CNode * m_prev;
  CNode * m_up;
  CNode * m_down;

  int m_level;
  int m_index;
  int m_stamp;
};


#endif
