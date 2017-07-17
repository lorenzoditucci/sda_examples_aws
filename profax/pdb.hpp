#ifndef _PDB_H
#define _PDB_H

#include <vector>
#include <fstream>
#include "MatVec.h"
#include <string.h>

using namespace std;

typedef vector<vector<REAL> > POSITIONS;

struct ATOM_ {
  char name[5];
  float pos[3];
};

struct AA {
  vector<ATOM_> atoms;
  int type;
};


void writeToPDB(const char * fname, const vector<AA> aas);
void writeLine(ofstream & fout, int index, const char * aname, const char * resname, 
	       const char * chainid, int resnum,
	       REAL x, REAL y, REAL z);
void loadCas(const char * fname, vector<vector<REAL> > & cas);


#endif
