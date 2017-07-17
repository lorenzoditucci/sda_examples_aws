#ifndef _RSS_H
#define _RSS_H

#include "bv.hpp"
#include <vector>
#include "MatVec.h"

class CRss : public CBV
{
public:
    // Construct an RSS to bound a list of vertices.
    CRss(const REAL vert[][3], int size)
    { computeRss(vert, size); }
    CRss(const CRss * rss1, const CRss * rss2,
         const REAL rot[3][3], const REAL trans[3]);
    
    
 //   virtual REAL computeDistance(const CBV * bv2,
   //                              const REAL rot[3][3], const REAL trans[3]);
    virtual void updateBV(const CBV * bv1, const CBV * bv2,
                          const REAL rot[3][3], const REAL trans[3]);
    
    virtual REAL getVolume()
    { return m_volume; }
    
    virtual void undoBV();
    
private:
    void computeRss(const REAL vert[][3], int size);
    void computeVertices();
    void extractVertices(const CRss * bv1, const CRss * bv2,
                         const REAL rot[3][3], const REAL trans[3],
                         REAL vert[16][3], int & size);
    
    // Compute the volume of the RSS.
    void computeVolume()
    { m_volume = m_radius*(2*m_sides[0]*m_sides[1] +
                           m_radius*M_PI*(m_sides[0]+m_sides[1]+
                                          (4.0/3.0)*m_radius)); }
    
    
    REAL m_vertices[8][3];
    
    REAL m_pos[3];       // position of rectangle
    REAL m_orient[3][3]; // orientation of rectangle
    REAL m_sides[2];     // length of its sides
    REAL m_radius;       // the radius of the swept sphere
    REAL m_volume;       // the volume of the RSS
    
    REAL m_undoPos[3];
    REAL m_undoOrient[3][3];
    REAL m_undoSides[2];
    REAL m_undoRadius;
    REAL m_undoVertices[8][3];
};


#endif

