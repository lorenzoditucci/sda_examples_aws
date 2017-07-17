#include <iomanip>
#include <fstream>
#include <cstdio>
#include <iostream>

#include "pdb.hpp"
#include "eef1.hpp"

void readline(const char * buf, char * aname, char * rname, int & resnum,
	      float & x, float & y, float & z);

// Load the coordinates of CA atoms from a PDB file.
void loadCas(const char * fname, POSITIONS & cas)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open input file " << fname << endl;
      exit(0);
    }

  char buf[200];
  char aname[5], rname[5];
  int resnum, curr_res = -1;;
  float x,y,z;
  while (!fin.eof())
    {
      fin.getline(buf, 199);
      if (strncmp(buf, "ATOM", 4) != 0)
	continue;

      readline(buf, aname, rname, resnum, x, y, z);
     
      if (strcmp("CA", aname) == 0)
	{
	  vector<REAL> a(3);
	  a[0] = x; a[1] = y; a[2] = z;  

	  cas.push_back(a);
	}
    }
}

// Read a line of atom coordinates from a PDB file.
void readline(const char * buf, char * aname, char * rname, int & resnum,
	      float & x, float & y, float & z)
{
  int i = 12;
  while (buf[i] == ' ')
    i++;

  int j = 0;
  while (buf[i] != ' ')
    aname[j++] = buf[i++];

  aname[j] = '\0';

  strncpy(rname, &(buf[17]),3);
  rname[3] = '\0';

  i = 22;
  while (buf[i] == ' ')
    i++;

  resnum = atoi(&(buf[i]));
  
  i = 30;
  while (buf[i] == ' ')
    i++;
  x = atof(&(buf[i]));

  i = 38;
  while (buf[i] == ' ')
    i++;
  y = atof(&(buf[i]));

  i = 46;
  while (buf[i] == ' ')
    i++;
  z = atof(&(buf[i]));
}

// Write a structure into a file in PDB format.
void writeToPDB(const char * fname, const vector<AA> aas)
{
  ofstream fout(fname);
  if (!fout.is_open())
    {
      cout << "Could not output coordinates to: " << fname << endl;
      return;

    }

  char buf[100];
  sprintf(buf,"HEADER    %-70s", fname);
  fout << buf << endl;
  sprintf(buf, "%-80s", "COMPND");
  fout << buf << endl;
  sprintf(buf, "%-80s", "SOURCE");
  fout << buf << endl;

  int i = 1, j = 1, n, m;
  vector<AA>::const_iterator aat = aas.begin();
  for (; aat != aas.end(); aat++)
    {
      int size = aat->atoms.size();
      assert(size >= 5);

      writeLine(fout, i++, aat->atoms[0].name, AA_NAMES[aat->type], " ", j, 
		aat->atoms[0].pos[0], aat->atoms[0].pos[1], aat->atoms[0].pos[2]);
      if (aat == aas.begin())
	writeLine(fout, i++, aat->atoms[4].name, AA_NAMES[aat->type], " ", j, 
		  aat->atoms[4].pos[0], aat->atoms[4].pos[1], aat->atoms[4].pos[2]);
      else
	writeLine(fout, i++, aat->atoms[2].name, AA_NAMES[aat->type], " ", j, 
		  aat->atoms[2].pos[0], aat->atoms[2].pos[1], aat->atoms[2].pos[2]);

      if (aat == aas.end()-1)
	{
	  writeLine(fout, i++, aat->atoms[size-3].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[size-3].pos[0], aat->atoms[size-3].pos[1], 
		    aat->atoms[size-3].pos[2]);
	  writeLine(fout, i++, aat->atoms[size-2].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[size-2].pos[0], aat->atoms[size-2].pos[1], 
		    aat->atoms[size-2].pos[2]);

	  m = size - 3;
	}
      else
	{
	  writeLine(fout, i++, aat->atoms[size-2].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[size-2].pos[0], aat->atoms[size-2].pos[1], 
		    aat->atoms[size-2].pos[2]);
	  writeLine(fout, i++, aat->atoms[size-1].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[size-1].pos[0], aat->atoms[size-1].pos[1], 
		    aat->atoms[size-1].pos[2]);

	  m = size - 2;
	}

      if (aat == aas.begin())
	{
	  writeLine(fout, i++, aat->atoms[1].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[1].pos[0], aat->atoms[1].pos[1], aat->atoms[1].pos[2]);
	  writeLine(fout, i++, aat->atoms[2].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[2].pos[0], aat->atoms[2].pos[1], aat->atoms[2].pos[2]);
	  writeLine(fout, i++, aat->atoms[3].name, AA_NAMES[aat->type], " ", j, 
		aat->atoms[3].pos[0], aat->atoms[3].pos[1], aat->atoms[3].pos[2]);

	  n = 5;
	}
      else
	{
	  writeLine(fout, i++, aat->atoms[1].name, AA_NAMES[aat->type], " ", j, 
		    aat->atoms[1].pos[0], aat->atoms[1].pos[1], aat->atoms[1].pos[2]);

	  n = 3;
	}

      for (int k = n; k < m; k++)
	writeLine(fout, i++, aat->atoms[k].name, AA_NAMES[aat->type], " ", j, 
		  aat->atoms[k].pos[0], aat->atoms[k].pos[1], aat->atoms[k].pos[2]);

      if (aat == aas.end()-1)
	writeLine(fout, i++, aat->atoms[size-1].name, AA_NAMES[aat->type], " ", j, 
		  aat->atoms[size-1].pos[0], aat->atoms[size-1].pos[1], 
		  aat->atoms[size-1].pos[2]);
	      
      j++;
    }

  sprintf(buf,"TER   %5d%6s%3s  %4d%54s",i,"",AA_NAMES[(aat-1)->type],j-1," ");
  fout << buf << endl;
  sprintf(buf, "MASTER    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%10s", 
	  0,0,0,0,0,0,0,0,i-1,1,0,0," ");  
  fout << buf << endl;
  sprintf(buf,"%-80s","END");
  fout << buf << endl;
}

// Write a line of one atoms coordinates into a PDB style file.
void writeLine(ofstream & fout, int index, const char * aname, const char * resname, 
	       const char * chainid, int resnum,
	       REAL x, REAL y, REAL z)
{
  fout << "ATOM  ";

  fout.width(5);
  //fout << right;
  fout.setf(ios::right, ios::adjustfield);
  fout << index;

  fout << ' ';

  fout.setf(ios::left, ios::adjustfield);
  if (aname[0] > '0' && aname[0] <= '9')
    {
      fout.width(4);
      fout << aname;
    }
  else
    {
      fout << ' ';
      fout.width(3);
      fout << aname;
    }

  fout << ' ';

  fout << resname;

  fout << ' ';

  fout << chainid;

  fout.setf(ios::right, ios::adjustfield);
  fout.width(4);
  fout << resnum;

  fout << "    ";
  
  char buf[50];
  sprintf(buf, "%8.3f%8.3f%8.3f%26s", x,y,z," ");
  fout << buf << endl;;
}
