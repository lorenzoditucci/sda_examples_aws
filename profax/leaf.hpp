#ifndef _LEAF_H
#define _LEAF_H

#include "eef1.hpp"
#include "bv.hpp"
#include "node.hpp"
#include <iostream>

// A structure describing one angle change.
struct ANGLE_CHANGE {
  int m_index;
  REAL m_angle;
};

struct ANGLE_CHANGE_COMP {
  bool operator() (const ANGLE_CHANGE & a, const ANGLE_CHANGE & b)
  {
    return (a.m_index < b.m_index);
  }
};

// A structure describing one rotamer change.
struct ROTAMER_CHANGE {
  int m_index;
  int m_rotIndex;
};

struct ROTAMER_CHANGE_COMP {
  bool operator() (const ROTAMER_CHANGE & a, const ROTAMER_CHANGE & b)
  {
    return (a.m_index < b.m_index);
  }
};

class CTerm;

class CLeaf : public CNode
{
public:
  // Create a leaf with the given rotamer value and type.
  CLeaf(ROTAMER::ROT_TYPE rType, int rotIndex, const REAL * joint, 
	int type, int index) : 	CNode(index), m_joint(joint), m_type(type), 
				m_rotType(rType), m_angle(0.0), m_energy(0.0), 
				m_torsionE(0.0) 
  { changeRotamer(rotIndex); }

  // Create a leaf given a set of coordinates for the atoms.
  CLeaf(int index, COORDS pos, int type) : CNode(index),
					   m_type(type), m_angle(0.0), m_energy(0.0),
					   m_torsionE(0.0), m_joint(NULL) {
	  //m_positions = pos;
	  memcpy(m_positions, pos, 12*3*sizeof(REAL));
  }
  
  void rotate(REAL angle);
  virtual void computePairEnergy(CNode * pNode, const REAL rot[3][3],
			      const REAL trans[3], CTerm * term,
			      bool bSeparated);
  virtual void computeSelfEnergy(CTerm * term);
  virtual bool findPairClash(CNode * pNode, const REAL rot[3][3],
			       const REAL trans[3], bool bSeparated);

  // A leaf is never in clash with itself.
  virtual bool findSelfClash()
  { return false; }
  virtual void undo();

  void changeRotamer(int rotIndex);

  // Accessors
  int getType() const
  { return m_type; }
  int getRotIndex() const
  { return m_rotIndex; }
  virtual bool isLeaf() const
  { return true; }
  CLeaf * getNext() const
  { return (CLeaf*) m_next; }
  
  COORDS  getPositions() const
  {
	  //std::cout << "m position is " << m_positions[0][0] << std::endl;
	  return m_positions;
  }

  REAL ** getPositionsReal() const
  {	return (REAL**)m_positions; }
  const REAL * getJoint() const
  { return m_joint; }
  REAL getAngle() const
  { return m_angle; }
  REAL getTorsionE() const
  { return m_torsionE; }
  REAL getTorsionChange() const
  { return m_torsionE - m_undoTorsionE; }

  static REAL m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];

private:
  void computeDistances(CLeaf * pLeaf, const REAL rot[3][3],
			const REAL trans[3]);

  const REAL * m_joint;
  REAL m_angle;

  int m_type;

  COORDS m_positions;
  int m_rotIndex;
  ROTAMER::ROT_TYPE m_rotType;
  int m_nPos;
  REAL m_energy;
  REAL m_torsionE;

  // undo variables
  REAL m_undoAngle;
  REAL m_undoTorsionE;
  int m_undoRotIndex;
};


#endif
