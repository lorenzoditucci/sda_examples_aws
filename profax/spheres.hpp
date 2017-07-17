#ifndef _SPHERES_H_
#define _SPHERES_H_

#include "bv.hpp"
#include <vector>
#include "MatVec.h"

class CSphere : public CBV
{
public:
    
    // Create the bounding sphere from a set of points.
    CSphere(const REAL vert[][3], int size)
    { computeSphere(vert, size); }
    CSphere(const CSphere * s1, const CSphere * s2,
            const REAL rot[3][3], const REAL trans[3])
    { computeSphere(s1, s2, rot, trans); }
    
    virtual void updateBV(const CBV * bv1, const CBV * bv2,
                          const REAL rot[3][3], const REAL trans[3]);
    
    virtual void undoBV();
    virtual REAL getVolume()
    { return 4.0/3.0*M_PI*m_rad*m_rad*m_rad; }
    
    // Compute the distance between two spheres.
    virtual REAL computeDistance(const CBV * bv2,
                                 const REAL rot[3][3], const REAL trans[3])
    {
        REAL T[3];
        REAL Ttemp[3];
        
        const CSphere * s1 = (CSphere*) this;
        const CSphere * s2 = (CSphere*) bv2;
        
        VmV(Ttemp, trans, s1->m_center);
        MxVpV(T, rot, s2->m_center, Ttemp);
        
        return Vlength(T) - (s1->m_rad + s2->m_rad);
    }
    
    REAL m_rad;
    REAL m_center[3];
    
    
private:
    void computeSphere(const REAL vert[][3], int size);
    void computeSphere(const CSphere * s1, const CSphere * s2,
                       const REAL rot[3][3], const REAL trans[3]);
    
    
    REAL m_undorad;
    REAL m_undocenter[3];		     
};


#endif
