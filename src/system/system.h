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

   System(string, vector<string>, string);

   int getNatoms() { return natoms; }
   const vector<AtomsRes>   getSystemData()  { return atoms; }

  ~System();


private:

   string gro_file;
   vector<string> itp_files;
   string top_file;

   bool water;

   int natoms;

   vector<AtomsRes> atoms;

   vector<int> resStart; 

   vector<string> molecules;
   vector<int> nmol;

   void readItp();
   void readGro();
   void readTop();
   void matchRes(const vector<AtomsRes>, vector<bool> &, int);
   void printAtomInfo(const vector<bool>);
   string exractAndTrim(const string &, const int, const int);

};

#endif
