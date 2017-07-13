#include "spheres.hpp"

// Compute a sphere to bound the given vertices
void CSphere::computeSphere(const REAL vert[][3], int size)
{
  // Only one vertex.
  if (size == 1)
    {
      m_rad = 0.0;
      VcV(m_center, vert[0]);
    }
  // Exactly two vertices
  else if (size == 2)
    {
      REAL v[3];
      VmV(v, vert[0], vert[1]);
      m_rad = Vlength(v)*0.5 + 0.001;

      VpV(v, vert[0], vert[1]);
      VxS(m_center, v, 0.5);
    }
  
  // More than two vertices. Simple but suboptimal method:
  // Center of sphere is mean of points, radius is distance from center to
  // furthest vertex.
  else
    {
      REAL cen[3];
      Videntity(cen);
   
      for (int i = 0; i < size; i++)
	{
	  cen[0] += vert[i][0];
	  cen[1] += vert[i][1];
	  cen[2] += vert[i][2];
	}

      VxS(m_center, cen, 1.0/size);

      REAL v[3];
      REAL max = -1;
      int ind = -1;
      for (int i = 0; i < size; i++)
	{
	  VmV(v, m_center, vert[i]);
	  REAL l = Vlength2(v);
	  if (l > max)
	    max = l;
	}

      m_rad = sqrt(max) + 0.001;
    }

  assert(m_rad >= 0.0);
}

// Compute a sphere to bound two given spheres.
void CSphere::computeSphere(const CSphere * s1, const CSphere * s2, 
			    const REAL rot[3][3], const REAL trans[3])
{
  REAL cen2[3];

  if (!s2)
    {
      VcV(m_center, s1->m_center);
      m_rad = s1->m_rad;
    }
  else
    {
      MxVpV(cen2, rot, s2->m_center, trans);
      assert(s1->m_rad >= 0.0 && s2->m_rad >= 0.0);

      REAL diff[3], sum[3], v1[3], v2[3];
      VmV(diff, cen2, s1->m_center);
      REAL l = Vlength(diff);

      // S2 is contained in S1
      if (l + s2->m_rad < s1->m_rad)
	{
	  VcV(m_center, s1->m_center);
	  m_rad = s1->m_rad;
	}

      // S1 is contained in S2
      else if (l + s1->m_rad < s2->m_rad)
	{
	  VcV(m_center, cen2);
	  m_rad = s2->m_rad;
	}

      // General case
      else
	{
	  m_rad = 0.5 * (l + s1->m_rad + s2->m_rad) + 0.001;
	  assert(m_rad >= 0.0);
  
	  VpV(sum, cen2, s1->m_center);
	  VxS(v1, diff, (s2->m_rad - s1->m_rad)/l);
	  VpV(v2, sum, v1);
	  VxS(m_center, v2, 0.5);
	}
      
      assert(m_rad >= s1->m_rad && m_rad >= s2->m_rad);
    }
}

// Update the sphere to bound the two given spheres.
void CSphere::updateBV(const CBV * bv1, const CBV * bv2,
		       const REAL rot[3][3], const REAL trans[3])
{ 
  m_undorad = m_rad;
  VcV(m_undocenter, m_center);
  
  computeSphere((CSphere*)bv1, (CSphere*)bv2, rot, trans); 
}

// Undo the changes to the sphere caused by the latest change.
void CSphere::undoBV()
{
  m_rad = m_undorad;
  VcV(m_center, m_undocenter);
}  
