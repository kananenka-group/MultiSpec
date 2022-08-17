#ifndef GRO_H
#define GRO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include "../util/util.h"
#include "../const/const.h"
#include "atomsRes.h"

using namespace std;

class System{

public:

   System(string, vector<string>);

   int getNatoms() { return natoms; }
   const vector<AtomsRes>   getSystemData()  { return atoms; }

  ~System();


private:

   string gro_file;
   vector<string> itp_files;

   bool water;

   int natoms;
   int nwater;
   int atoms_in_water;

   vector<AtomsRes> atoms;

   waterM W;

   // all residue numbers (natoms)
   //vector<int> aResN;
   // unique residue numbers
   //vector<int> uResN;
   // all residue names (natoms)
   //vector<string> aRes;
   // unique residue names
   //vector<string> uRes;
   // all atom name (natoms)
   //vector<string> aAtoms;
   // unique atoms
   //vector<string> uAtoms;
   // atoms that belong to each residue name
   //vector<vector<int>> resAtoms;
   // starting indices for each residue
   vector<int> resStart; 
   vector<int> wOind;

   void readItp();
   void readGro();
   void readAtoms();
   void common();
   string exractAndTrim(const string &, const int, const int);

};

#endif
