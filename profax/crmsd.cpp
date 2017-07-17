#include "crmsd.hpp"
#include "pdb.hpp"

#include <iostream>

#define X 0
#define Y 1
#define Z 2
#define W 3


// Compute the cRMSD between the structures stored in the two given files.
REAL computeRMSD(const char * file1, const char * file2, REAL rot[3][3], REAL trans[3])
{
  vector<vector<REAL> > cas1, cas2;
  loadCas(file1, cas1);
  loadCas(file2, cas2);

  if (cas1.size() != cas2.size())
    {
      cout << "The two compared structures are of different length!!!" << endl;
      exit(0);
    }

  REAL res = CRMSD(cas1, cas2, rot, trans);

  return res;
}

REAL CRMSD(const vector<vector<REAL> > &  A, const vector<vector<REAL> > & B,
	   REAL rot[3][3], REAL trans[3])
{
  REAL Bt[A.size()][3];

  compute_alignment(A, B, rot, trans);

  // Transform the second chain to optimally align with the first.
  for (int k = 0; k < A.size(); k++)
    {
      Bt[k][X] = B[k][X] * rot[0][0] + B[k][Y] * rot[1][0] +
	B[k][Z] * rot[2][0] + trans[0];
      Bt[k][Y] = B[k][X] * rot[0][1] + B[k][Y] * rot[1][1] +
	B[k][Z] * rot[2][1] + trans[1];
      Bt[k][Z] = B[k][X] * rot[0][2] + B[k][Y] * rot[1][2] +
	B[k][Z] * rot[2][2] + trans[2];
    }

  REAL rmsd = 0;
  for (int i = 0; i < A.size(); i++) 
    {
      REAL a0 = Bt[i][X] - A[i][X];
      REAL a1 = Bt[i][Y] - A[i][Y];
      REAL a2 = Bt[i][Z] - A[i][Z];
      
      rmsd += (a0*a0 + a1*a1 + a2*a2); 
    }

  return sqrt(rmsd / A.size());
}

void compute_alignment(const vector<vector<REAL> > &  A, const vector<vector<REAL> > & B,
		       REAL rot[3][3], REAL trans[3])
{
  int i,j,k;
  REAL c1[3],c2[3];   /* center of mass for two point collections */
  REAL v1[3],v2[3];
  REAL recip;
  REAL tr;
  REAL m[4][4], q[4][4];
  REAL v[4];
  REAL cov[3][3];
  REAL aij[3][3];
  REAL quat[4];

  /* find the center of mass for the two collections */
  c1[X] = c1[Y] = c1[Z] = 0;
  c2[X] = c2[Y] = c2[Z] = 0;

  for (int i = 0; i < A.size(); i++) 
    {
      c1[X] += A[i][X];
      c1[Y] += A[i][Y];
      c1[Z] += A[i][Z];
      
      c2[X] += B[i][X];
      c2[Y] += B[i][Y];
      c2[Z] += B[i][Z];
    }
  
  recip = 1.0 / A.size();

  c1[X] *= recip;
  c1[Y] *= recip;
  c1[Z] *= recip;

  c2[X] *= recip;
  c2[Y] *= recip;
  c2[Z] *= recip;

  /* create the cross-covariance matrix */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      cov[i][j] = 0;

  for (int i = 0; i < A.size(); i++)
    {
      v1[X] = A[i][X] - c1[X];
      v1[Y] = A[i][Y] - c1[Y];
      v1[Z] = A[i][Z] - c1[Z];
      
      v2[X] = B[i][X] - c2[X];
      v2[Y] = B[i][Y] - c2[Y];
      v2[Z] = B[i][Z] - c2[Z];
      
      cov[X][X] += v1[X] * v2[X];
      cov[X][Y] += v1[X] * v2[Y];
      cov[X][Z] += v1[X] * v2[Z];
      
      cov[Y][X] += v1[Y] * v2[X];
      cov[Y][Y] += v1[Y] * v2[Y];
      cov[Y][Z] += v1[Y] * v2[Z];
      
      cov[Z][X] += v1[Z] * v2[X];
      cov[Z][Y] += v1[Z] * v2[Y];
      cov[Z][Z] += v1[Z] * v2[Z];
    }

  /* aij = cov - transpose(cov) */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      aij[i][j] = cov[i][j] - cov[j][i];

  /* find the trace of the covariance matrix */
  tr = cov[X][X] + cov[Y][Y] + cov[Z][Z];

  m[0][0] = tr;

  m[1][0] = m[0][1] = aij[1][2];
  m[2][0] = m[0][2] = aij[2][0];
  m[3][0] = m[0][3] = aij[0][1];

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m[i+1][j+1] = cov[i][j] + cov[j][i] - (i == j) * tr;

  /* find the eigenvector corresponding to the largest eigenvalue */
  /* of this matrix */
  Meigen4(q, v, m);

  if( v[0] > v[1] ){
      if( v[0] > v[2] ){
	  if( v[0] > v[3] )
	    i = 0;
	  else
	    i = 3;
      }
      else{
	  if( v[2] > v[3] )
	    i = 2;
	  else
	    i =3;
      }
  }
  else{
      if( v[1] > v[2] ){
	  if( v[1] > v[3] )
	    i = 1;
	  else
	    i = 3;
      }
      else{
	  if( v[2] > v[3] )
	    i = 2;
	  else
	    i =3;
      }
  }
  
  quat[0] = q[0][i];
  quat[1] = q[1][i];
  quat[2] = q[2][i];
  quat[3] = q[3][i];
  
  convert_quat_to_mat(quat, rot);

  REAL c3[3];

  /* determine best translation */
  vapply (rot, c2, c3);
  trans[0] = c1[0] - c3[0];
  trans[1] = c1[1] - c3[1];
  trans[2] = c1[2] - c3[2];
}


void convert_quat_to_mat(REAL q[4], REAL mat[3][3])
{
 double q00,q01,q02,q03;
 double q11,q12,q13;
 double q22,q23;
 double q33;

 q00 = q[0] * q[0];
 q01 = q[0] * q[1];
 q02 = q[0] * q[2];
 q03 = q[0] * q[3];

 q11 = q[1] * q[1];
 q12 = q[1] * q[2];
 q13 = q[1] * q[3];

 q22 = q[2] * q[2];
 q23 = q[2] * q[3];

 q33 = q[3] * q[3];

 mat[X][X] = q00 + q11 - q22 - q33;
 mat[X][Y] = 2 * (q12 - q03);
 mat[X][Z] = 2 * (q13 + q02);

 mat[Y][X] =  2 * (q12 + q03);
 mat[Y][Y] =  q00 + q22 - q11 - q33;
 mat[Y][Z] = 2 * (q23 - q01);

 mat[Z][X] = 2 * (q13 - q02);
 mat[Z][Y] = 2 * (q23 + q01);
 mat[Z][Z] = q00 + q33 - q11 - q22;
}


void vapply (REAL m[3][3], const REAL a[3], REAL b[3])
{
  int j;

  for (j = 0; j <= 2; j++)
    b[j] = a[0] * m[0][j] + a[1] * m [1][j] +
      a[2] * m[2][j];
}

