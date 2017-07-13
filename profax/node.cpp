#include <cmath>
#include <cassert>
#include <cfloat>

#include "node.hpp"
#include "slist.hpp"
#include "pairtree.hpp"
#include "rss.hpp"
#include "spheres.hpp"
#include "node.hpp"

extern REAL CUTOFF_DISTANCE;
vector<CTerm*> CNode::m_undoPairs;

// Create a new internal node above the node that is given
CNode::CNode(CNode * down, CSkiplist * pSL)
: m_next(NULL), m_prev(NULL), m_up(NULL), m_index(down->getIndex()),
m_down(down), m_bv(NULL), m_stamp(0), m_level(down->getLevel()+1),
m_nc(BOTH), m_pSL(pSL)
{
    Midentity(m_rotate);
    down->setUp(this);
    CBV* next_bv = down->getNext() ? down->getNext()->m_bv : NULL;
    
    // Create a bounding volume for the new node.
#ifdef USE_RSS
    m_bv = (CBV*) new CRss((CRss*)(down->m_bv), (CRss*)(next_bv),
                           down->m_rotate, down->m_translate);
#else
    m_bv = (CBV*) new CSphere((CSphere*)(down->m_bv),
                              (CSphere*)(next_bv),
                              down->m_rotate, down->m_translate);
#endif
}

// Compute the transform to the next node at the same level.
// The transform is the product of the transforms held by the two
// child nodes of this node.
void CNode::transformToNext()
{
    assert(!isLeaf());
    assert(goDown());
    
    VcV(m_undoTranslate, m_translate);
    McM(m_undoRotate, m_rotate);
    
    CNode * pNode = goDown();
    if (!pNode->getNext())
        return;
    
    MxVpV(m_translate, pNode->m_rotate, pNode->getNext()->m_translate,
          pNode->m_translate);
    MxM(m_rotate, pNode->m_rotate, pNode->getNext()->m_rotate);
}

// Compute the energy between this node and the given node.
// We are also given the transform between the corrdinate frames
// of the two nodes.
void CNode::computePairEnergy(CNode * pNode, const REAL rot[3][3],
                              const REAL trans[3], CTerm * term,
                              bool bSeparated)
{
    assert(this != pNode);
    int time = m_pSL->getTime();;
    
    // Do not perform distance comp. for the top 3 levels (a heuristic)
    REAL dist = -1.0;
    if (getLevel() < m_pSL->getLevels() - 3)
        dist = m_bv->computeDistance(pNode->getBV(), rot, trans);
    
    // Proceed down the hierarchy only if distance is smaller than the
    // cutoff.
    if (dist <= CUTOFF_DISTANCE)
    {
        CNode * son1_1 = goDown();
        CNode * son1_2 = son1_1->getNext();
        CNode * son2_1 = pNode->goDown();
        CNode * son2_2 = son2_1->getNext();
        
        REAL rot_[3][3];
        REAL trans_[3];
        REAL Ttemp[3];
        REAL Rtemp[3][3];
        bool bSep;
        bool bTest;
        
        // Decide whether the first two sons are separated.
        bSep = bSeparated || son1_2->isSeparator(time);
        
        // The first sons should be tested if they are separated or either is
        // affected.
        bTest = bSep || son1_1->isSeparator(time) || son1_1->isAffected(time) ||
        son2_1->isAffected(time);
        if (bTest)
            son1_1->computePairEnergy(son2_1, rot, trans,
                                      term->getChild(0), bSep);
        
        // update the rotation and translation.
        VmV(Ttemp, trans, son1_1->m_translate);
        MTxV(trans_, son1_1->m_rotate, Ttemp);
        
        MTxM(rot_, son1_1->m_rotate, rot);
        
        // Decide whether the second son of this node and the first of te other
        // should be tested.
        bTest = bSeparated || son1_2->isSeparator(time) || son1_2->isAffected(time) ||
        son2_1->isAffected(time);
        if (bTest)
            son1_2->computePairEnergy(son2_1, rot_, trans_,
                                      term->getChild(1), bSeparated);
        
        // The second son of teh other node may not exist if we are at the
        // right end of the hierarchy.
        if (son2_2)
        {
            MxVpV(trans_, rot, son2_1->m_translate, trans);
            MxM(rot_, rot, son2_1->m_rotate);
            
            bSep = bSeparated  || son1_2->isSeparator(time) ||
            son2_1->isSeparator(time);
            bTest = bSep || son1_1->isSeparator(time) || son1_1->isAffected(time) ||
            son2_2->isAffected(time);
            if (bTest)
                son1_1->computePairEnergy(son2_2, rot_, trans_,
                                          term->getChild(2), bSep);
            
            // update the rotation and translation.
            VmV(Ttemp, trans_, son1_1->m_translate);
            MTxV(trans_, son1_1->m_rotate, Ttemp);
            
            McM(Rtemp, rot_);
            MTxM(rot_, son1_1->m_rotate, Rtemp);
            
            bSep = bSeparated || son2_1->isSeparator(time);
            bTest = bSep || son1_2->isSeparator(time) || son1_2->isAffected(time) ||
            son2_2->isAffected(time);
            if (bTest)
                son1_2->computePairEnergy(son2_2, rot_, trans_,
                                          term->getChild(3), bSep);
        }
        
        // Recompute the interaction energy of the two nodes since it changed.
        term->recompute(son2_2 != NULL);
    }
    else
    {
        // The two nodes are too far to interact reset their interaction energy.
        term->reset();
    }
    
    // Save a pointer to this node in case we need to undo the last move.
    m_undoPairs.push_back(term);
    
    return;
}

// Compute the energy contributed by interactions inside this node.
void CNode::computeSelfEnergy(CTerm * term)
{
    assert(!isLeaf());
    
    int time = m_pSL->getTime();
    CNode * curr = goDown();
    
    // Compute the energy inside the first child node.
    if (curr->isAffected(time))
        curr->computeSelfEnergy(term->getChild(0));
    
    CNode * next = curr->getNext();
    if (next)
    {
        // Compute the energy inside the second child.
        curr->computePairEnergy(next, curr->m_rotate,
                                curr->m_translate,
                                term->getChild(1),
                                false);
        
        // Test the second child for self collisions.
        if (next->isAffected(time))
            next->computeSelfEnergy(term->getChild(2));
    }
    
    // Recompute the energy inside this node.
    term->recompute();
    
    // Save a pointer to this node in case we need to undo the last move.
    m_undoPairs.push_back(term);
    
    return;
}


bool CNode::findPairClash(CNode * pNode, const REAL rot[3][3],
                          const REAL trans[3], bool bSeparated)
{
    assert(!isLeaf());
    assert(this != pNode);
    int time = m_pSL->getTime();
    
    REAL dist = m_bv->computeDistance(pNode->getBV(), rot, trans);
    
    if (dist < CLASH_CUTOFF_DISTANCE)
    {
        CNode * son1_1 = goDown();
        CNode * son1_2 = son1_1->getNext();
        CNode * son2_1 = pNode->goDown();
        CNode * son2_2 = son2_1->getNext();
        
        REAL rot_[3][3];
        REAL trans_[3];
        REAL Ttemp[3];
        REAL Rtemp[3][3];
        bool bSep;
        bool bTest;
        
        // The first sons are separated if the second son
        // of the left node is marked.
        bSep = bSeparated || (son1_2->isSeparator(time));
        bTest = bSep || son1_1->isAffected(time) ||  son1_1->isSeparator(time) ||
        son2_1->isAffected(time);
        if (bTest && son1_1->findPairClash(son2_1, rot, trans, bSep))
            return true;
        
        // update the rotation and translation.
        VmV(Ttemp, trans, son1_1->m_translate);
        MTxV(trans_, son1_1->m_rotate, Ttemp);
        
        MTxM(rot_, son1_1->m_rotate, rot);
        
        bTest = bSeparated || son1_2->isSeparator(time) || son1_2->isAffected(time) ||
        son2_1->isAffected(time);
        if (bTest && son1_2->findPairClash(son2_1, rot_, trans_, bSeparated))
            return true;
        
        if (son2_2)
        {
            MxVpV(trans_, rot, son2_1->m_translate, trans);
            MxM(rot_, rot, son2_1->m_rotate);
            
            bSep = bSeparated  || son1_2->isSeparator(time) ||
            son2_1->isSeparator(time);
            bTest = bSep || son1_1->isSeparator(time) || son1_1->isAffected(time) ||
            son2_2->isAffected(time);
            if (bTest && son1_1->findPairClash(son2_2, rot_, trans_, bSep))
                return true;
            
            // update the rotation and translation.
            VmV(Ttemp, trans_, son1_1->m_translate);
            MTxV(trans_, son1_1->m_rotate, Ttemp);
            
            McM(Rtemp, rot_);
            MTxM(rot_, son1_1->m_rotate, Rtemp);
            
            bSep = bSeparated || son2_1->isSeparator(time);
            bTest = bSep || son1_2->isSeparator(time) || son1_2->isAffected(time) ||
            son2_2->isAffected(time);
            if (bTest && son1_2->findPairClash(son2_2, rot_, trans_, bSep))
                return true;
        }
    }
    
    return false;
}

bool CNode::findSelfClash()
{
    assert(!isLeaf());
    
    int time = m_pSL->getTime();
    CNode * curr = goDown();
    
    // Find clashes inside the first child.
    if (curr->isAffected(time) && curr->findSelfClash())
        return true;
    
    CNode * next = curr->getNext();
    if (next)
    {
        // Find clashes between the two children.
        if (curr->findPairClash(next, curr->m_rotate, 
                                curr->m_translate, 
                                false))
            return true;
        
        // Find clashes inside the second child.
        if (next->isAffected(time) && next->findSelfClash())
            return true;
    }
    
    return false;
}

// Recompute the bounding volume to bound the two BVs of the child nodes 
// below this one in the hierarchy.
void CNode::updateBV()
{
    assert(!isLeaf());
    
    CNode * down = goDown();
    CBV* next_bv = down->getNext() ? down->getNext()->m_bv : NULL;
    
    m_bv->updateBV(down->m_bv, next_bv, 
                   down->m_rotate, down->m_translate);
}

// Undo any changes made to this node during the last move.
void CNode::undo()
{
    assert(!isLeaf());
    
    if (m_nc & CNode::TRANSFORM)
    {
        McM(m_rotate, m_undoRotate);
        VcV(m_translate, m_undoTranslate);
    }
    
    if (m_nc & CNode::BOX)
    {
        undoBV();
    }
}






