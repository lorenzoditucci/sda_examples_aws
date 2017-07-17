#include <vector>
#include "rss.hpp"
#include "RectDist.h"

// Compute the 8 vertices of the bounding box of this RSS
void CRss::computeVertices() 
{
  REAL rcol1[3], rcol2[3], rcol3[3];
  McolcV(rcol1, m_orient, 0);
  McolcV(rcol2, m_orient, 1);
  McolcV(rcol3, m_orient, 2);

  REAL rtemp[3] = {m_radius, m_radius, m_radius};
  REAL temp[3];

  // Compute the first vertex
  MxV(temp, m_orient, rtemp);
  VmV(m_vertices[0], m_pos, temp);
 
  // compute the second vertex
  VpVxS(m_vertices[1], m_vertices[0], rcol3, 2*m_radius);

  // compute the third vertex
  VpVxS(m_vertices[2], m_vertices[0], rcol1, 2*m_radius + m_sides[0]);

  // compute the fourth vertex
  VpVxS(m_vertices[3], m_vertices[2], rcol3, 2*m_radius);

  // compute the fifth vertex
  VpVxS(m_vertices[4], m_vertices[0], rcol2, 2*m_radius + m_sides[1]);

  // compute the sixth vertex
  VpVxS(m_vertices[5], m_vertices[4], rcol3, 2*m_radius);

  // compute the seventh vertex
  VpVxS(m_vertices[6], m_vertices[4], rcol1, 2*m_radius + m_sides[0]);

  // compute the eighth vertex
  VpVxS(m_vertices[7], m_vertices[6], rcol3, 2*m_radius);
}

// Get the vertices of the two given RSSs
void CRss::extractVertices(const CRss * rss1, const CRss * rss2,
			   const REAL rot[3][3], const REAL trans[3],
			   REAL vert[16][3], int & size)
{
  // The first RSS has the same reference frame as me.
  for (int i = 0; i < 8; i++)
    VcV(vert[i], rss1->m_vertices[i]);

  if (rss2)
    {
      // Add the vertices to the list.
      for (int k = 0; k < 8; k++) 
	MxVpV(vert[8+k], rot, rss2->m_vertices[k], trans);
	
      size = 16;
    } 
  else
    size = 8;
}

// Construct an RSS to bound the two given RSSs.
CRss::CRss(const CRss * rss1, const CRss * rss2, 
	   const REAL rot[3][3], const REAL trans[3])
{
  REAL vertices[16][3];
  int size;
  
  extractVertices(rss1, rss2, rot, trans, vertices, size);
  computeRss(vertices, size);
}

// Update the RSS to bound the two given RSSs
void CRss::updateBV(const CBV * bv1, const CBV * bv2,
		    const REAL rot[3][3], const REAL trans[3])
{
  REAL vertices[16][3];
  int size;

  CRss * rss1 = (CRss*) bv1;
  CRss * rss2 = (CRss*) bv2;

  extractVertices(rss1, rss2, rot, trans, vertices, size);

  // Save the current RSS in case we need to undo
  VcV(m_undoPos, m_pos);
  m_undoRadius = m_radius;
  McM(m_undoOrient, m_orient);
  m_undoSides[0] = m_sides[0];
  m_undoSides[1] = m_sides[1];
  
  for (int i = 0; i < 8; i++)
    VcV(m_undoVertices[i], m_vertices[i]);  

  computeRss(vertices, size);
}

#define MaxOfTwo(a,b) ((a)>(b) ? (a) : (b))

// Compute an RSS to bound a list of vertices (taken from the PQP library).
void CRss::computeRss(const REAL vert[][3], int size)
{
  for (int i = 0; i < 8; i++)
    VcV(m_undoVertices[i], m_vertices[i]);  

  REAL cov[3][3];
  computeCov(vert, size, cov);
  
  // Compute the eigenvectors and sort in decreasing order.
  REAL eigval[3], eigvec[3][3];
  Meigen(eigvec, eigval, cov);
  sortEV(m_orient, eigvec, eigval);
  
  // Project all points onto new coordinates
  REAL points[size][3];
  for (int i = 0; i < size; i++)
    MTxV(points[i], m_orient, vert[i]);
  
  REAL minx, maxx, miny, maxy, minz, maxz, c[3];

  // compute thickness, which determines radius, and z of 
  // rectangle corner
  REAL cz,radsqr;
  minz = maxz = points[0][2];
  for (int i = 1; i < size; i++) 
    {
      if (points[i][2] < minz) 
	minz = points[i][2];
      else if (points[i][2] > maxz) 
	maxz = points[i][2];
    }

  m_radius = (REAL)0.5*(maxz - minz);
  
  // If the points are coplanar set the radius to 0 (itay)
  if (m_radius < 1e-10)
    m_radius = 0.0;
  
  radsqr = m_radius*m_radius;
  cz = (REAL)0.5*(maxz + minz);

  // compute an initial length of rectangle along x direction

  // find minx and maxx as starting points
  int minindex, maxindex;
  minindex = maxindex = 0;
  for (int i = 1; i < size; i++) 
  {
    if (points[i][0] < points[minindex][0]) 
      minindex = i; 
    else if (points[i][0] > points[maxindex][0]) 
      maxindex = i;
  }

  REAL x, dz;
  dz = points[minindex][2] - cz;
  minx = points[minindex][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = points[maxindex][2] - cz;
  maxx = points[maxindex][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow minx
  for (int i = 0; i < size; i++) 
  {
    if (points[i][0] < minx) 
    {
      dz = points[i][2] - cz;
      x = points[i][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x < minx) 
	minx = x;
    }
  }

  // grow maxx
  for (int i = 0; i < size; i++) 
  {
    if (points[i][0] > maxx) 
    {
      dz = points[i][2] - cz;
      x = points[i][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x > maxx) 
	maxx = x;
    }
  }
  
  // compute an initial length of rectangle along y direction
  // find miny and maxy as starting points
  minindex = maxindex = 0;
  for (int i = 1; i < size; i++) 
  {
    if (points[i][1] < points[minindex][1]) 
      minindex = i;
    else if (points[i][1] > points[maxindex][1]) 
      maxindex = i;
  }
  
  REAL y;
  dz = points[minindex][2] - cz;
  miny = points[minindex][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = points[maxindex][2] - cz;
  maxy = points[maxindex][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow miny
  for (int i = 0; i < size; i++) 
    {
      if (points[i][1] < miny) 
	{
	  dz = points[i][2] - cz;
	  y = points[i][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
	  if (y < miny) 
	    miny = y;
	}
    }
  
  // grow maxy
  for (int i = 0; i < size; i++) 
    {
      if (points[i][1] > maxy) 
	{
	  dz = points[i][2] - cz;
	  y = points[i][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
	  if (y > maxy) 
	    maxy = y;
	}
    }
  
  // corners may have some points which are not covered - 
  // grow lengths if necessary
  REAL dx, dy, u, t;
  REAL a = sqrt((REAL)0.5);
  for (int i = 0; i < size; i++) 
    {
      if (points[i][0] > maxx) 
	{
	  if (points[i][1] > maxy) 
	    {
	      dx = points[i][0] - maxx;
	      dy = points[i][1] - maxy;
	      u = dx*a + dy*a;
	      t = (a*u - dx)*(a*u - dx) + 
		(a*u - dy)*(a*u - dy) +
		(cz - points[i][2])*(cz - points[i][2]);
	      u = u - sqrt(MaxOfTwo(radsqr - t,0));
	      if (u > 0) 
		{
		  maxx += u*a;
		  maxy += u*a;
		}
	    }
	  else if (points[i][1] < miny) 
	    {
	      dx = points[i][0] - maxx;
	      dy = points[i][1] - miny;
	      u = dx*a - dy*a;
	      t = (a*u - dx)*(a*u - dx) + 
		(-a*u - dy)*(-a*u - dy) +
		(cz - points[i][2])*(cz - points[i][2]);
	      u = u - sqrt(MaxOfTwo(radsqr - t,0));
	      if (u > 0) 
		{
		  maxx += u*a;
		  miny -= u*a;
		}
	    }
	}
      else if (points[i][0] < minx) 
	{
	  if (points[i][1] > maxy) 
	    {
	      dx = points[i][0] - minx;
	      dy = points[i][1] - maxy;
	      u = dy*a - dx*a;
	      t = (-a*u - dx)*(-a*u - dx) + 
		(a*u - dy)*(a*u - dy) +
		(cz - points[i][2])*(cz - points[i][2]);
	      u = u - sqrt(MaxOfTwo(radsqr - t,0));
	      if (u > 0) 
		{
		  minx -= u*a;
		  maxy += u*a;
		}     
	    }
	  else if (points[i][1] < miny) 
	    {
	      dx = points[i][0] - minx;
	      dy = points[i][1] - miny;
	      u = -dx*a - dy*a;
	      t = (-a*u - dx)*(-a*u - dx) + 
		(-a*u - dy)*(-a*u - dy) +
		(cz - points[i][2])*(cz - points[i][2]);
	      u = u - sqrt(MaxOfTwo(radsqr - t,0));
	      if (u > 0) 
		{
		  minx -= u*a; 
		  miny -= u*a;
		}
	    }
	}
    }

  c[0] = minx;
  c[1] = miny;
  c[2] = cz;
  MxV(m_pos, m_orient, c);
  
  // Save the length of both sides of the rectangle
  m_sides[0] = maxx - minx;  
  if (m_sides[0] < 0) 
    m_sides[0] = 0;
  
  m_sides[1] = maxy - miny;
  if (m_sides[1] < 0)
    m_sides[1] = 0;
  
  computeVertices();
  computeVolume();
}

// Undo the changes to this RSS casued by the latest move.
void CRss::undoBV()
{
  VcV(m_pos, m_undoPos);
  m_radius = m_undoRadius;
  McM(m_orient, m_undoOrient);
  m_sides[0] = m_undoSides[0];
  m_sides[1] = m_undoSides[1];

  for (int i = 0; i < 8; i++)
    VcV(m_vertices[i], m_undoVertices[i]);

  computeVolume();
}

// Compute the minimal distance between this RSS and the given RSS.
/*REAL CRss::computeDistance(const CBV * bv2,
			   const REAL rot[3][3], const REAL trans[3])
{
  REAL T[3];
  REAL R[3][3];

  CRss * rss1 = (CRss*) this;
  CRss * rss2 = (CRss*) bv2;

  REAL Ttemp[3];
  VmV(Ttemp, trans, m_pos);
  MxVpV(Ttemp, rot, rss2->m_pos, Ttemp);
  MTxV(T, rss1->m_orient, Ttemp);
 
  REAL Rtemp[3][3];
  MxM(Rtemp, rot, rss2->m_orient);
  MTxM(R, rss1->m_orient, Rtemp);

  REAL dist = sqrt(RectDist(R, T, rss1->m_sides, rss2->m_sides)) - 
    (rss1->m_radius + rss2->m_radius);

  return dist;
}*/



