#ifndef _CRMSD_H_
#define _CRMSD_H_

#include <vector>

#include "MatVec.h"

#include <stdlib.h>

using namespace std;

void vapply (REAL m[3][3], const REAL a[3], REAL b[3]);
void convert_quat_to_mat(REAL q[4], REAL mat[3][3]);
void compute_alignment(const vector<vector<REAL> > & A, const vector<vector<REAL> > & B,
		       REAL rot[3][3], REAL trans[3]);

REAL CRMSD(const vector<vector<REAL> > & A, const vector<vector<REAL> > & B,
	   REAL rot[3][3], REAL trans[3]);

REAL computeRMSD(const char * file1, const char * file2, REAL rot[3][3], REAL trans[3]);

#endif
