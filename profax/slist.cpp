#include <cmath>
#include <cassert>
#include <cfloat>
#include <list>

#include "slist.hpp"

// Create a hierarchy of bounding volumes over the chain.
void CSkiplist::createHierarchy()
{
  assert(m_chain.getLength() > 0);

  CNode * currLevel = (CNode*) m_chain.getLink(0);
  CNode * nextLevel;
  CNode * currNode, *prevNode = NULL, *pNode = currLevel;

  m_levels = 1;

  // Create new levels as long as the current level has more 
  // than one node
  while(currLevel->getNext() != NULL)
    {
      while(pNode != NULL)
	{
	  currNode = new CNode(pNode, this);
	  assert(currNode);

	  // Set the transformation to the next node
	  currNode->transformToNext();

	  // Save the pointer to the first node of the new level being created
	  if (pNode == currLevel)
	    nextLevel = currNode;

	  // Adjust pointers
	  currNode->setPrev(prevNode);
	  if (prevNode)
	    prevNode->setNext(currNode);

	  prevNode = currNode;

	  // Advance the pointer on the lower level two steps forward 
	  // if possible.
	  pNode = pNode->getNext();
	  if (pNode)
	    pNode = pNode->getNext();
	} 

      currLevel = nextLevel;
      prevNode = NULL;
      pNode = currLevel;

      m_levels++;
    }
  
  m_root = currLevel;

  // Create the energy tree structure
  m_pairTree = new CTerm(m_levels-1, CTerm::IDENTICAL);

  // Set the pointer to the hierarchy for each link in the chain.
  for (int i = 0; i < getLength(); i++)
    getLink(i)->m_pSL = this;
}

// Apply a move to the chaintree.
// A move consists of changes to backbone angles and to rotamers.
void CSkiplist::makeMove(const vector<ANGLE_CHANGE> & angles,
			 const vector<ROTAMER_CHANGE> & rotamers)
{
  // Clear the list of nodes to undo.
  clearUndo();

  CNode* nodes[angles.size() + rotamers.size()];
  CNode::NODE_CHANGE updates[angles.size() + rotamers.size()];

  // Increment the timer to mark the beginning of a new time step.
  incTime();

  int a = 0, b = 0, i = 0, k = 0;
  // Mark the affected leaf nodes according to the type of change they undergo.
  // The changes are handled in order from left to right.
  while(angles[a].m_index < getLength() || rotamers[b].m_index < getLength())
    { 
      CLeaf * pLeaf;

      // The next change is only an angle change.
      if (angles[a].m_index < rotamers[b].m_index)
	{
	  pLeaf = getLink(angles[a].m_index);
	  updates[i] = CNode::TRANSFORM;

	  pLeaf->rotate(angles[a].m_angle);

	  saveUndoNode((CNode*) pLeaf);
	  a++;
	}

      // The next change is only a rotamer change.
      else if (angles[a].m_index > rotamers[b].m_index)
      {
	pLeaf = getLink(rotamers[b].m_index);
	updates[i] = CNode::BOX;

	pLeaf->changeRotamer(rotamers[b].m_rotIndex);

	saveUndoNode(pLeaf);
	b++;
      }
      
      // The next change is both an angle and a rotamer change.
      else
      {
	assert(angles[a].m_index == rotamers[b].m_index);

	pLeaf = getLink(angles[a].m_index);
	updates[i] = CNode::BOTH;

	pLeaf->rotate(angles[a].m_angle);	
	pLeaf->changeRotamer(rotamers[b].m_rotIndex);
	
	saveUndoNode(pLeaf);
	a++; b++;
      }

      // Set the time stamp on the affected leaf node
      pLeaf->setStamp(m_time, updates[i]);
	
      // Go up to the parent.
      CNode * parent = pLeaf->getParent();
    
      // If two nodes converge (have the same parent) then save
      // only the first one.
      if (k == 0 || parent != nodes[k-1])
	{
	  nodes[k] = parent;
	  updates[k] = updates[i];

	  if (pLeaf->goUp())
	    updates[k] = (CNode::NODE_CHANGE) (updates[k] | CNode::BOX);

	  k++;
	}
      else
	updates[k-1] = (CNode::NODE_CHANGE) (updates[k-1] | updates[i]);

      i++;
    }

  // Update the rest of the levels of the hierarchy.
  int num = k;
  while (nodes[0] != m_root)
    {
      k = 0;

      for (i = 0; i < num; i++)
	{
	  CNode * pNode = nodes[i];
	  
	  // The node is marked to signal that it has changed.
	  pNode->setStamp(m_time, updates[i]);
	  saveUndoNode(pNode);
	  
	  // Update the transformation to the next node.
	  if ((updates[i] & CNode::TRANSFORM) &&  pNode->getNext())
	    pNode->transformToNext();
	 
	  // Recompute the BV if necessary.
   	  if (updates[i] & CNode::BOX)
	    pNode->updateBV();

	  // Go up to the parent
	  CNode * parent = pNode->getParent();

	  // Update the parent only if it is not the parent 
	  // of the previous node (converging paths)
	  if (k == 0 || parent != nodes[k-1])
	    {
	      nodes[k] = parent;
	      updates[k] = updates[i];
	      
	      if (pNode->goUp())
		updates[k] = (CNode::NODE_CHANGE) (updates[k] | CNode::BOX);

	      k++;
	    }
	  else
	    updates[k-1] = (CNode::NODE_CHANGE) (updates[k-1] | updates[i]);
	}
	
	num = k;
    }

  // The root box was not updated. There is no need to actually recompute the 
  // BV of the root, only update its time stamp.
  assert(nodes[0] == m_root);
  nodes[0]->setStamp(m_time, CNode::BOX);   
}

// Clear the undo lists of the chaintree and the energytree.
void CSkiplist::clearUndo()
{
  m_undoNodes.clear();
  assert(m_undoNodes.empty());

  CNode::m_undoPairs.clear();
  assert(CNode::m_undoPairs.empty());
}

// Undo the last move that was applied to the data structures.
// Undos both chaintree and energytree
void CSkiplist::undoLastMove()
{
  vector<CNode*>::iterator nit = m_undoNodes.begin();
  for (; nit != m_undoNodes.end(); nit++)
    (*nit)->undo();

  vector<CTerm*>::iterator pit = CNode::m_undoPairs.begin();
  for (; pit != CNode::m_undoPairs.end(); pit++)
    (*pit)->undo();

  // Also undo torsion energy computation.
  m_chain.undoTorsionChange();
  m_energy = m_undoEnergy;
}

REAL CSkiplist::computeEnergy(const vector<ANGLE_CHANGE> & angles)
{
  m_root->computeSelfEnergy(m_pairTree);

  m_undoEnergy = m_energy;

  m_energy = m_pairTree->m_energy + m_chain.updateTorsionEnergy(angles) +
    m_chain.getSolventE(); 

  return m_energy;
}

// Find a steric clash in the protein. Used a filter before complete
// energy computation.
bool CSkiplist::findSelfClash()
{
  return m_root->findSelfClash();
}

// An auxiliary function to compute the current distance between two atoms. 
// The  atoms are specified by the indices of their respective links and the 
// index of each atom in its respective link.
REAL CSkiplist::computeDistance(int ind1, int i, int ind2, int j)
{
  CLeaf * pLeaf1, * pLeaf2;
  pLeaf1 = getLink(ind1);
  pLeaf2 = getLink(ind2);

  // Get the position of each atom in its leaf's coordinate frame.
  REAL trans1[3], trans2[3], temp[3];
  VcV(trans1, pLeaf1->getPositions()[i]);
  VcV(trans2, pLeaf2->getPositions()[j]);

  CNode * pNode1, * pNode2;

  pNode1 = (CNode*) pLeaf1;
  pNode2 = (CNode*) pLeaf2;

  // Recurtsively climb up the tree until the paths from both leaves converge.
  // Update the translation from the current node to both leaf as you climb.
  while (pNode1 != pNode2)
    {
      // Go up a level on the first path.
      if (!pNode1->goUp())
	{
	  assert(pNode1->getPrev());
	  pNode1 = pNode1->getPrev();
	  MxVpV(temp, pNode1->m_rotate, trans1, pNode1->m_translate);	  
	  VcV(trans1, temp);
	}
     
      assert(pNode1->goUp());
      pNode1 = pNode1->goUp();

      // Go up a level on the second path.
      if (!pNode2->goUp())
	{
	  assert(pNode2->getPrev());
	  pNode2 = pNode2->getPrev();
	  MxVpV(temp, pNode2->m_rotate, trans2, pNode2->m_translate);
	  VcV(trans2, temp);
	}
     
      assert(pNode2->goUp());
      pNode2 = pNode2->goUp();
    }

  // Compute the difference between the two translations and return its norm.
  VmV(temp, trans1, trans2);
  return Vlength(temp);
}

