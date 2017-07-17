/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             S. Gottschalk
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#ifndef PQP_MATVEC_H
#define PQP_MATVEC_H


#include <cmath>
#include <cstdio>
#include <cassert>

typedef float REAL;
typedef REAL (*COORDS)[3];
//typedef REAL COORDS[12][3];


#ifndef M_PI
const REAL M_PI = (REAL)3.14159265359;
#endif

#ifdef gnu
#include "zzzz.h"

#ifdef hppa
#define myfabs(x) \
 ({double __value, __arg = (x); \
   asm("fabs,dbl %1, %0": "=f" (__value): "f" (__arg)); \
   __value; \
});
#endif

#ifdef mips
#define myfabs(x) \
 ({double __value, __arg = (x); \
   asm("abs.d %0, %1": "=f" (__value): "f" (__arg)); \
   __value; \
});
#endif

#else  

#define myfabs(x) ((x < 0) ? -x : x)

#endif

//inline float sqrt(float x) { return (float)sqrt((double)x); }
//inline float cos(float x) { return (float)cos((double)x); }
//inline float sin(float x) { return (float)sin((double)x); }
//inline float fabs(float x) { return (float)fabs((double)x); }

inline
void
Mprintg(const REAL M[3][3])
{
  printf("%g %g %g\n%g %g %g\n%g %g %g\n",
	 M[0][0], M[0][1], M[0][2],
	 M[1][0], M[1][1], M[1][2],
	 M[2][0], M[2][1], M[2][2]);
}


inline
void
Mfprint(FILE *f, const REAL M[3][3])
{
  fprintf(f, "%g %g %g\n%g %g %g\n%g %g %g\n",
	 M[0][0], M[0][1], M[0][2],
	 M[1][0], M[1][1], M[1][2],
	 M[2][0], M[2][1], M[2][2]);
}

inline
void
Mprint(const REAL M[3][3])
{
  printf("%g %g %g\n%g %g %g\n%g %g %g\n",
	 M[0][0], M[0][1], M[0][2],
	 M[1][0], M[1][1], M[1][2],
	 M[2][0], M[2][1], M[2][2]);
}

inline
void
Mprint4(const REAL M[4][4])
{
  printf("%g %g %g % g\n%g %g %g %g\n%g %g %g %g\n%g %g %g %g\n",
	 M[0][0], M[0][1], M[0][2], M[0][3],
	 M[1][0], M[1][1], M[1][2], M[1][3],
	 M[2][0], M[2][1], M[2][2], M[2][3],
	 M[3][0], M[3][1], M[3][2], M[3][3]);
}

inline
void
Vprintg(const REAL V[3])
{
  printf("%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Vfprint(FILE *f, const REAL V[3])
{
  fprintf(f, "%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Vprint(const REAL V[3])
{
  printf("%g %g %g\n", V[0], V[1], V[2]);
}

inline
void
Vprint4(const REAL V[4])
{
  printf("%g %g %g %g\n", V[0], V[1], V[2], V[3]);
}

inline
void
Midentity(REAL M[3][3])
{
  M[0][0] = M[1][1] = M[2][2] = 1.0;
  M[0][1] = M[1][2] = M[2][0] = 0.0;
  M[0][2] = M[1][0] = M[2][1] = 0.0;
}

inline
void
Midentity4(REAL M[4][4])
{
  M[0][0] = M[1][1] = M[2][2] = M[3][3] = 1.0;
  M[0][1] = M[0][2] = M[0][3] = 0.0;
  M[1][0] = M[1][2] = M[1][3] = 0.0;
  M[2][0] = M[2][1] = M[2][3] = 0.0;
  M[3][0] = M[3][1] = M[3][2] = 0.0;
}

inline
void
Videntity(REAL T[3])
{
  T[0] = T[1] = T[2] = 0.0;
}

inline
void
McM(REAL Mr[3][3], const REAL M[3][3])
{
  Mr[0][0] = M[0][0];  Mr[0][1] = M[0][1];  Mr[0][2] = M[0][2];
  Mr[1][0] = M[1][0];  Mr[1][1] = M[1][1];  Mr[1][2] = M[1][2];
  Mr[2][0] = M[2][0];  Mr[2][1] = M[2][1];  Mr[2][2] = M[2][2];
}

inline
void
MTcM(REAL Mr[3][3], const REAL M[3][3])
{
  Mr[0][0] = M[0][0];  Mr[1][0] = M[0][1];  Mr[2][0] = M[0][2];
  Mr[0][1] = M[1][0];  Mr[1][1] = M[1][1];  Mr[2][1] = M[1][2];
  Mr[0][2] = M[2][0];  Mr[1][2] = M[2][1];  Mr[2][2] = M[2][2];
}

inline
void
VcV(REAL Vr[3], const REAL V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

inline
void
McolcV(REAL Vr[3], const REAL M[3][3], int c)
{
  Vr[0] = M[0][c];
  Vr[1] = M[1][c];
  Vr[2] = M[2][c];
}

inline
void
VcolcM(REAL Mr[3][3], const REAL V[3], int c)
{
  Mr[0][c] = V[0];
  Mr[1][c] = V[1];
  Mr[2][c] = V[2];
}

inline
void
McolcMcol(REAL Mr[3][3], int cr, const REAL M[3][3], int c)
{
  Mr[0][cr] = M[0][c];
  Mr[1][cr] = M[1][c];
  Mr[2][cr] = M[2][c];
}

inline
void
MxMpV(REAL Mr[3][3], const REAL M1[3][3], const REAL M2[3][3], const REAL T[3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
	      M1[0][1] * M2[1][0] +
	      M1[0][2] * M2[2][0] +
	      T[0]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
	      M1[1][1] * M2[1][0] +
	      M1[1][2] * M2[2][0] +
	      T[1]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
	      M1[2][1] * M2[1][0] +
	      M1[2][2] * M2[2][0] +
	      T[2]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
	      M1[0][1] * M2[1][1] +
	      M1[0][2] * M2[2][1] +
	      T[0]);
  Mr[1][1] = (M1[1][0] * M2[0][1] +
	      M1[1][1] * M2[1][1] +
 	      M1[1][2] * M2[2][1] +
	      T[1]);
  Mr[2][1] = (M1[2][0] * M2[0][1] +
	      M1[2][1] * M2[1][1] +
	      M1[2][2] * M2[2][1] +
	      T[2]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
	      M1[0][1] * M2[1][2] +
	      M1[0][2] * M2[2][2] +
	      T[0]);
  Mr[1][2] = (M1[1][0] * M2[0][2] +
	      M1[1][1] * M2[1][2] +
	      M1[1][2] * M2[2][2] +
	      T[1]);
  Mr[2][2] = (M1[2][0] * M2[0][2] +
	      M1[2][1] * M2[1][2] +
	      M1[2][2] * M2[2][2] +
	      T[2]);
}

inline
void
MxM(REAL Mr[3][3], const REAL M1[3][3], const REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
	      M1[0][1] * M2[1][0] +
	      M1[0][2] * M2[2][0]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
	      M1[1][1] * M2[1][0] +
	      M1[1][2] * M2[2][0]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
	      M1[2][1] * M2[1][0] +
	      M1[2][2] * M2[2][0]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
	      M1[0][1] * M2[1][1] +
	      M1[0][2] * M2[2][1]);
  Mr[1][1] = (M1[1][0] * M2[0][1] +
	      M1[1][1] * M2[1][1] +
 	      M1[1][2] * M2[2][1]);
  Mr[2][1] = (M1[2][0] * M2[0][1] +
	      M1[2][1] * M2[1][1] +
	      M1[2][2] * M2[2][1]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
	      M1[0][1] * M2[1][2] +
	      M1[0][2] * M2[2][2]);
  Mr[1][2] = (M1[1][0] * M2[0][2] +
	      M1[1][1] * M2[1][2] +
	      M1[1][2] * M2[2][2]);
  Mr[2][2] = (M1[2][0] * M2[0][2] +
	      M1[2][1] * M2[1][2] +
	      M1[2][2] * M2[2][2]);
}


inline
void
MxMT(REAL Mr[3][3], const REAL M1[3][3], const REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
	      M1[0][1] * M2[0][1] +
	      M1[0][2] * M2[0][2]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
	      M1[1][1] * M2[0][1] +
	      M1[1][2] * M2[0][2]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
	      M1[2][1] * M2[0][1] +
	      M1[2][2] * M2[0][2]);
  Mr[0][1] = (M1[0][0] * M2[1][0] +
	      M1[0][1] * M2[1][1] +
	      M1[0][2] * M2[1][2]);
  Mr[1][1] = (M1[1][0] * M2[1][0] +
	      M1[1][1] * M2[1][1] +
 	      M1[1][2] * M2[1][2]);
  Mr[2][1] = (M1[2][0] * M2[1][0] +
	      M1[2][1] * M2[1][1] +
	      M1[2][2] * M2[1][2]);
  Mr[0][2] = (M1[0][0] * M2[2][0] +
	      M1[0][1] * M2[2][1] +
	      M1[0][2] * M2[2][2]);
  Mr[1][2] = (M1[1][0] * M2[2][0] +
	      M1[1][1] * M2[2][1] +
	      M1[1][2] * M2[2][2]);
  Mr[2][2] = (M1[2][0] * M2[2][0] +
	      M1[2][1] * M2[2][1] +
	      M1[2][2] * M2[2][2]);
}

inline
void
MTxM(REAL Mr[3][3], const REAL M1[3][3], const REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
	      M1[1][0] * M2[1][0] +
	      M1[2][0] * M2[2][0]);
  Mr[1][0] = (M1[0][1] * M2[0][0] +
	      M1[1][1] * M2[1][0] +
	      M1[2][1] * M2[2][0]);
  Mr[2][0] = (M1[0][2] * M2[0][0] +
	      M1[1][2] * M2[1][0] +
	      M1[2][2] * M2[2][0]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
	      M1[1][0] * M2[1][1] +
	      M1[2][0] * M2[2][1]);
  Mr[1][1] = (M1[0][1] * M2[0][1] +
	      M1[1][1] * M2[1][1] +
 	      M1[2][1] * M2[2][1]);
  Mr[2][1] = (M1[0][2] * M2[0][1] +
	      M1[1][2] * M2[1][1] +
	      M1[2][2] * M2[2][1]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
	      M1[1][0] * M2[1][2] +
	      M1[2][0] * M2[2][2]);
  Mr[1][2] = (M1[0][1] * M2[0][2] +
	      M1[1][1] * M2[1][2] +
	      M1[2][1] * M2[2][2]);
  Mr[2][2] = (M1[0][2] * M2[0][2] +
	      M1[1][2] * M2[1][2] +
	      M1[2][2] * M2[2][2]);
}

inline
void
MxV(REAL Vr[3], const REAL M1[3][3], const REAL V1[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
	   M1[0][1] * V1[1] + 
	   M1[0][2] * V1[2]);
  Vr[1] = (M1[1][0] * V1[0] +
	   M1[1][1] * V1[1] + 
	   M1[1][2] * V1[2]);
  Vr[2] = (M1[2][0] * V1[0] +
	   M1[2][1] * V1[1] + 
	   M1[2][2] * V1[2]);
}


inline
void
MxVpV(REAL Vr[3], const REAL M1[3][3], const REAL V1[3], const REAL V2[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
	   M1[0][1] * V1[1] + 
	   M1[0][2] * V1[2] + 
	   V2[0]);
  Vr[1] = (M1[1][0] * V1[0] +
	   M1[1][1] * V1[1] + 
	   M1[1][2] * V1[2] + 
	   V2[1]);
  Vr[2] = (M1[2][0] * V1[0] +
	   M1[2][1] * V1[1] + 
	   M1[2][2] * V1[2] + 
	   V2[2]);
}


inline
void
sMxVpV(REAL Vr[3], REAL s1, const REAL M1[3][3], const REAL V1[3], const REAL V2[3])
{
  Vr[0] = s1 * (M1[0][0] * V1[0] +
		M1[0][1] * V1[1] + 
		M1[0][2] * V1[2]) +
		V2[0];
  Vr[1] = s1 * (M1[1][0] * V1[0] +
		M1[1][1] * V1[1] + 
		M1[1][2] * V1[2]) + 
		V2[1];
  Vr[2] = s1 * (M1[2][0] * V1[0] +
		M1[2][1] * V1[1] + 
		M1[2][2] * V1[2]) + 
		V2[2];
}

inline
void
MTxV(REAL Vr[3], const REAL M1[3][3], const REAL V1[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
	   M1[1][0] * V1[1] + 
	   M1[2][0] * V1[2]); 
  Vr[1] = (M1[0][1] * V1[0] +
	   M1[1][1] * V1[1] + 
	   M1[2][1] * V1[2]);
  Vr[2] = (M1[0][2] * V1[0] +
	   M1[1][2] * V1[1] + 
	   M1[2][2] * V1[2]); 
}

inline
void
sMTxV(REAL Vr[3], REAL s1, const REAL M1[3][3], const REAL V1[3])
{
  Vr[0] = s1*(M1[0][0] * V1[0] +
	      M1[1][0] * V1[1] + 
	      M1[2][0] * V1[2]); 
  Vr[1] = s1*(M1[0][1] * V1[0] +
	      M1[1][1] * V1[1] + 
	      M1[2][1] * V1[2]);
  Vr[2] = s1*(M1[0][2] * V1[0] +
	      M1[1][2] * V1[1] + 
	      M1[2][2] * V1[2]); 
}

inline
void
sMxV(REAL Vr[3], REAL s1, const REAL M1[3][3], const REAL V1[3])
{
  Vr[0] = s1*(M1[0][0] * V1[0] +
	      M1[0][1] * V1[1] + 
	      M1[0][2] * V1[2]); 
  Vr[1] = s1*(M1[1][0] * V1[0] +
	      M1[1][1] * V1[1] + 
	      M1[1][2] * V1[2]);
  Vr[2] = s1*(M1[2][0] * V1[0] +
	      M1[2][1] * V1[1] + 
	      M1[2][2] * V1[2]); 
}


inline
void
VmV(REAL Vr[3], const REAL V1[3], const REAL V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

inline
void
VpV(REAL Vr[3], const REAL V1[3], const REAL V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

inline
void
VpVxS(REAL Vr[3], const REAL V1[3], const REAL V2[3], REAL s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

inline
void
VmVxS(REAL Vr[3], const REAL V1[3], const REAL V2[3], REAL s)
{
  Vr[0] = V1[0] - V2[0] * s;
  Vr[1] = V1[1] - V2[1] * s;
  Vr[2] = V1[2] - V2[2] * s;
}

inline 
void
MskewV(REAL M[3][3], const REAL v[3])
{
  M[0][0] = M[1][1] = M[2][2] = 0.0;
  M[1][0] = v[2];
  M[0][1] = -v[2];
  M[0][2] = v[1];
  M[2][0] = -v[1];
  M[1][2] = -v[0];
  M[2][1] = v[0];
}


inline
void
VcrossV(REAL Vr[3], const REAL V1[3], const REAL V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

inline
REAL
Vlength(const REAL V[3])
{
  return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

inline
REAL
Vlength2(const REAL V[3])
{
  return (V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

inline
REAL 
Vdist(const REAL V1[3], const REAL V2[3])
{
  return sqrt((V1[0]-V2[0])*(V1[0]-V2[0]) + 
	      (V1[1]-V2[1])*(V1[1]-V2[1]) + 
	      (V1[2]-V2[2])*(V1[2]-V2[2]));
}

inline
REAL 
Vdist2(const REAL V1[3], const REAL V2[3])
{
  return ((V1[0]-V2[0])*(V1[0]-V2[0]) + 
	  (V1[1]-V2[1])*(V1[1]-V2[1]) + 
	  (V1[2]-V2[2])*(V1[2]-V2[2]));
}

inline
void
Vnormalize(REAL V[3])
{
  REAL d = (REAL)1.0 / sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
  V[0] *= d;
  V[1] *= d;
  V[2] *= d;
}

inline
REAL
VdotV(const REAL V1[3], const REAL V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

inline
REAL
VdistV2(const REAL V1[3], const REAL V2[3])
{
  return ( (V1[0]-V2[0]) * (V1[0]-V2[0]) + 
	   (V1[1]-V2[1]) * (V1[1]-V2[1]) + 
	   (V1[2]-V2[2]) * (V1[2]-V2[2]));
}

inline 
void
MpVout(REAL Mr[3][3], const REAL V1[3])
{
  Mr[0][0] += V1[0] * V1[0]; 
  Mr[0][1] += V1[0] * V1[1];
  Mr[0][2] += V1[0] * V1[2];
  Mr[1][1] += V1[1] * V1[1];
  Mr[1][2] += V1[1] * V1[2];
  Mr[2][2] += V1[2] * V1[2];
}



inline
void
VxS(REAL Vr[3], const REAL V[3], REAL s)
{
  Vr[0] = V[0] * s;
  Vr[1] = V[1] * s;
  Vr[2] = V[2] * s;
}

inline
void
MRotZ(REAL Mr[3][3], REAL t)
{
  Mr[0][0] = cos(t);
  Mr[1][0] = sin(t);
  Mr[0][1] = -Mr[1][0];
  Mr[1][1] = Mr[0][0];
  Mr[2][0] = Mr[2][1] = 0.0;
  Mr[0][2] = Mr[1][2] = 0.0;
  Mr[2][2] = 1.0;
}

inline
void
MRotX(REAL Mr[3][3], REAL t)
{
  Mr[1][1] = cos(t);
  Mr[2][1] = sin(t);
  Mr[1][2] = -Mr[2][1];
  Mr[2][2] = Mr[1][1];
  Mr[0][1] = Mr[0][2] = 0.0;
  Mr[1][0] = Mr[2][0] = 0.0;
  Mr[0][0] = 1.0;
}

inline
void
MRotY(REAL Mr[3][3], REAL t)
{
  Mr[2][2] = cos(t);
  Mr[0][2] = sin(t);
  Mr[2][0] = -Mr[0][2];
  Mr[0][0] = Mr[2][2];
  Mr[1][2] = Mr[1][0] = 0.0;
  Mr[2][1] = Mr[0][1] = 0.0;
  Mr[1][1] = 1.0;
}

inline
void
MVtoOGL(double oglm[16], const REAL R[3][3], const REAL T[3])
{
  oglm[0] = (double)R[0][0]; 
  oglm[1] = (double)R[1][0]; 
  oglm[2] = (double)R[2][0]; 
  oglm[3] = 0.0;
  oglm[4] = (double)R[0][1]; 
  oglm[5] = (double)R[1][1];
  oglm[6] = (double)R[2][1];
  oglm[7] = 0.0;
  oglm[8] = (double)R[0][2];
  oglm[9] = (double)R[1][2];
  oglm[10] = (double)R[2][2];
  oglm[11] = 0.0;
  oglm[12] = (double)T[0];
  oglm[13] = (double)T[1];
  oglm[14] = (double)T[2];
  oglm[15] = 1.0;
}

inline 
void
OGLtoMV(REAL R[3][3], REAL T[3], const double oglm[16])
{
  R[0][0] = (REAL)oglm[0];
  R[1][0] = (REAL)oglm[1];
  R[2][0] = (REAL)oglm[2];

  R[0][1] = (REAL)oglm[4];
  R[1][1] = (REAL)oglm[5];
  R[2][1] = (REAL)oglm[6];

  R[0][2] = (REAL)oglm[8];
  R[1][2] = (REAL)oglm[9];
  R[2][2] = (REAL)oglm[10];

  T[0] = (REAL)oglm[12];
  T[1] = (REAL)oglm[13];
  T[2] = (REAL)oglm[14];
}

// taken from quatlib, written by Richard Holloway
const int QX = 0;
const int QY = 1;
const int QZ = 2;
const int QW = 3;

inline
void 
MRotQ(REAL destMatrix[3][3], REAL srcQuat[4])
{
  REAL  s;
  REAL  xs, ys, zs,
    	    wx, wy, wz,
	        xx, xy, xz,
	        yy, yz, zz;

  /* 
   * For unit srcQuat, just set s = 2.0; or set xs = srcQuat[QX] + 
   *   srcQuat[QX], etc. 
   */

  s = (REAL)2.0 / (srcQuat[QX]*srcQuat[QX] + srcQuat[QY]*srcQuat[QY] + 
    	     srcQuat[QZ]*srcQuat[QZ] + srcQuat[QW]*srcQuat[QW]);

  xs = srcQuat[QX] * s;   ys = srcQuat[QY] * s;   zs = srcQuat[QZ] * s;
  wx = srcQuat[QW] * xs;  wy = srcQuat[QW] * ys;  wz = srcQuat[QW] * zs;
  xx = srcQuat[QX] * xs;  xy = srcQuat[QX] * ys;  xz = srcQuat[QX] * zs;
  yy = srcQuat[QY] * ys;  yz = srcQuat[QY] * zs;  zz = srcQuat[QZ] * zs;

  destMatrix[QX][QX] = (REAL)1.0 - (yy + zz);
  destMatrix[QX][QY] = xy + wz;
  destMatrix[QX][QZ] = xz - wy;

  destMatrix[QY][QX] = xy - wz;
  destMatrix[QY][QY] = (REAL)1.0 - (xx + zz);
  destMatrix[QY][QZ] = yz + wx;

  destMatrix[QZ][QX] = xz + wy;
  destMatrix[QZ][QY] = yz - wx;
  destMatrix[QZ][QZ] = (REAL)1.0 - (xx + yy);
} 

inline
void
Mqinverse(REAL Mr[3][3], REAL m[3][3])
{
  int i,j;

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
    {
      int i1 = (i+1)%3;
      int i2 = (i+2)%3;
      int j1 = (j+1)%3;
      int j2 = (j+2)%3;
      Mr[i][j] = (m[j1][i1]*m[j2][i2] - m[j1][i2]*m[j2][i1]);
    }
}

// Meigen from Numerical Recipes in C
#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

inline
void
Meigen(REAL vout[3][3], REAL dout[3], REAL a[3][3])
{
  int n = 3;
  int j,iq,ip,i;
  REAL tresh,theta,tau,t,sm,s,h,g,c;
  int nrot;
  REAL b[3];
  REAL z[3];
  REAL v[3][3];
  REAL d[3];
  
  Midentity(v);
  for(ip=0; ip<n; ip++) 
    {
      b[ip] = a[ip][ip];
      d[ip] = a[ip][ip];
      z[ip] = 0.0;
    }
  
  nrot = 0;
  
  for(i=0; i<50; i++)
    {

      sm=0.0;
      for(ip=0;ip<n;ip++) for(iq=ip+1;iq<n;iq++) sm+=fabs(a[ip][iq]);
      if (sm == 0.0)
	{
	  McM(vout, v);
	  VcV(dout, d);
	  return;
	}
      
      
      if (i < 3) tresh=(REAL)0.2*sm/(n*n);
      else tresh=0.0;
      
      for(ip=0; ip<n; ip++) for(iq=ip+1; iq<n; iq++)
	{
	  g = (REAL)100.0*fabs(a[ip][iq]);
	  if (i>3 && 
	      fabs(d[ip])+g==fabs(d[ip]) && 
	      fabs(d[iq])+g==fabs(d[iq]))
	    a[ip][iq]=0.0;
	  else if (fabs(a[ip][iq])>tresh)
	    {
	      h = d[iq]-d[ip];
	      if (fabs(h)+g == fabs(h)) t=(a[ip][iq])/h;
	      else
		{
		  theta=(REAL)0.5*h/(a[ip][iq]);
		  t=(REAL)(1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
		  if (theta < 0.0) t = -t;
		}
	      c=(REAL)1.0/sqrt(1+t*t);
	      s=t*c;
	      tau=s/((REAL)1.0+c);
	      h=t*a[ip][iq];
	      z[ip] -= h;
	      z[iq] += h;
	      d[ip] -= h;
	      d[iq] += h;
	      a[ip][iq]=0.0;
	      for(j=0;j<ip;j++) { ROTATE(a,j,ip,j,iq); } 
	      for(j=ip+1;j<iq;j++) { ROTATE(a,ip,j,j,iq); } 
	      for(j=iq+1;j<n;j++) { ROTATE(a,ip,j,iq,j); } 
	      for(j=0;j<n;j++) { ROTATE(v,j,ip,j,iq); } 
	      nrot++;
	    }
	}
      for(ip=0;ip<n;ip++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }

  fprintf(stderr, "eigen: too many iterations in Jacobi transform.\n");

  return;
}

inline
void
Meigen4(REAL v[4][4], REAL d[4], REAL a[4][4])
{
  int n = 4;
  int j,iq,ip,i;
  REAL tresh,theta,tau,t,sm,s,h,g,c;
  int nrot;
  REAL b[4];
  REAL z[4];
  
  Midentity4(v);
  for(ip=0; ip<n; ip++) 
    {
      b[ip] = a[ip][ip];
      d[ip] = a[ip][ip];
      z[ip] = 0.0;
    }
  
  nrot = 0;
  
  for(i=0; i<50; i++)
    {

      sm=0.0;
      for(ip=0;ip<n;ip++) for(iq=ip+1;iq<n;iq++) sm+=fabs(a[ip][iq]);
      if (sm == 0.0)
	  return;
      
      if (i < 3) 
	tresh=(REAL)0.2*sm/(n*n);
      else 
	tresh=0.0;
      
      for(ip=0; ip<=n-2; ip++) for(iq=ip+1; iq<=n-1; iq++)
	{
	  g = (REAL)100.0*fabs(a[ip][iq]);
	  if (i>3 && 
	      fabs(d[ip])+g==fabs(d[ip]) && 
	      fabs(d[iq])+g==fabs(d[iq]))
	    a[ip][iq]=0.0;
	  else if (fabs(a[ip][iq])>tresh)
	    {
	      h = d[iq]-d[ip];
	      if (fabs(h)+g == fabs(h)) t=(a[ip][iq])/h;
	      else
		{
		  theta=(REAL)0.5*h/(a[ip][iq]);
		  t=(REAL)(1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
		  if (theta < 0.0) t = -t;
		}
	      c=(REAL)1.0/sqrt(1+t*t);
	      s=t*c;
	      tau=s/((REAL)1.0+c);
	      h=t*a[ip][iq];
	      z[ip] -= h;
	      z[iq] += h;
	      d[ip] -= h;
	      d[iq] += h;
	      a[ip][iq]=0.0;
	      for(j=0;j<ip;j++) { ROTATE(a,j,ip,j,iq); } 
	      for(j=ip+1;j<iq;j++) { ROTATE(a,ip,j,j,iq); } 
	      for(j=iq+1;j<n;j++) { ROTATE(a,ip,j,iq,j); } 
	      for(j=0;j<n;j++) { ROTATE(v,j,ip,j,iq); } 
	      nrot++;
	    }
	}
      for(ip=0;ip<=n-1;ip++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }

  fprintf(stderr, "eigen: too many iterations in Jacobi transform.\n");

  return;
}

inline
void
sortEV(REAL out[3][3], const REAL in[3][3], const REAL val[3])
{
  // place axes of in matrix in order of increasing eigen values
  int min, mid, max;

  if (val[0] > val[1]) 
    { 
      max = 0; min = 1; 
    }
  else 
    {
      min = 0; max = 1; 
    }
    
  if (val[2] < val[min]) 
    { 
      mid = min; min = 2; 
    }
  else if (val[2] > val[max])
    {
      mid = max; max = 2; 
    }
  else 
    { 
      mid = 2; 
    }
    
  assert(val[max] >= val[mid] && val[max] >= val[min] && val[mid] >= val[min]);

  McolcMcol(out, 0, in, max);
  McolcMcol(out, 1, in, mid);

  out[0][2] = in[1][max]*in[2][mid] - in[1][mid]*in[2][max];
  out[1][2] = in[0][mid]*in[2][max] - in[0][max]*in[2][mid];
  out[2][2] = in[0][max]*in[1][mid] - in[0][mid]*in[1][max];
}

inline
void
compute_rotation(REAL R[3][3], const REAL axis[3], REAL angle)
{
  REAL kxkx = axis[0] * axis[0];
  REAL kyky = axis[1] * axis[1];
  REAL kzkz = axis[2] * axis[2];

  REAL kxky = axis[0] * axis[1];
  REAL kxkz = axis[0] * axis[2];
  REAL kykz = axis[1] * axis[2];

  REAL c = cos(angle);
  REAL s = sin(angle);
  REAL v = 1 - c;

  R[0][0] = kxkx * v + c;
  R[1][1] = kyky * v + c;
  R[2][2] = kzkz * v + c;

  R[0][1] = kxky * v - axis[2]*s;
  R[1][0] = kxky * v + axis[2]*s;

  R[0][2] = kxkz * v + axis[1]*s;
  R[2][0] = kxkz * v - axis[1]*s;
  
  R[1][2] = kykz * v - axis[0]*s;
  R[2][1] = kykz * v + axis[0]*s;
}

inline
void
computeCov(const REAL vert[][3], int size, REAL cov[3][3])
{
  cov[0][0] = cov[1][1] = cov[2][2] = 0.0;
  cov[0][1] = cov[1][2] = cov[2][0] = 0.0;
  cov[0][2] = cov[1][0] = cov[2][1] = 0.0;

  REAL mean[3];
  Videntity(mean);
  for (int i = 0; i < size; i++)
    VpV(mean, mean, vert[i]);

  VxS(mean, mean, 1/((double) size));

  REAL temp[3]; 
  for (int i = 0; i < size; i++)
    {
      VmV(temp, vert[i], mean);
      MpVout(cov, temp);
    }

  cov[1][0] = cov[0][1];
  cov[2][1] = cov[1][2];
  cov[2][0] = cov[0][2];
}

inline
REAL
compute_angle(const REAL x[3], const REAL y[3], const REAL z[3]) 
{
  REAL	acc, d1, d2;
  
  acc = 0;
  for (int i = 0; i < 3; i++)
    acc += (y[i] - x[i]) * (y[i] - z[i]);
  
  d1 = Vdist(x, y);
  d2 = Vdist(z, y);
  
  acc = acc / (d1 * d2);
  
  if (acc > 1.0)
    acc = 1.0;
  else if (acc < -1.0)
    acc = -1.0;
  
  acc = acos(acc);
  
  return acc;
}

inline
REAL
compute_dihed(const REAL a[3], const REAL b[3], const REAL c[3], const REAL d[3])
{
  REAL	q[3], r[3], s[3], t[3], u[3], v[3], z[3];
  REAL acc, w;
  
  Videntity(z);

  for (int i = 0; i < 3; i++) {
    q[i] = b[i] - a[i];
    r[i] = b[i] - c[i];
    s[i] = c[i] - d[i];
  }
  
  VcrossV(t, q, r);
  VcrossV(u, s, r);
  VcrossV(v, u, t);
  
  w = VdotV(v, r);
  
  acc = compute_angle(t, z, u);
  
  if (w < 0) 
    acc = -acc;

  return (acc);
}


#endif
// MATVEC_H
