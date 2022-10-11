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
   const vector<int> getChgSt()  { return chgSt; }
   const vector<int> getCgInd()  { return grpInd; }

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
   vector<int> chgSt;
   vector<int> ChgAll;
   vector<int> grpInd;

   void readItp();
   void readGro();
   void readTop();
   void matchRes(const vector<AtomsRes>, vector<bool> &, int, const vector<int>,
                 vector<int> &);
   void printAtomInfo(const vector<bool>);
   string exractAndTrim(const string &, const int, const int);

};

#endif
