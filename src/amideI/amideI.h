#ifndef AMIDE_H
#define AMIDE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm> // for transform
#include <iterator> // begin, end, ostream_iterator
#include <numeric>  // iota
#include "../system/system.h"
#include "../traj/traj.h"
#include "../vec/vecManip.h"
#include "../const/const.h"
#include "../util/util.h"

using namespace std;

class amideI{

public:
   amideI(string, string, vector<string>, string, string,
          vector<string>);

  ~amideI();

private:

  System s;

  string gro_file;
  string traj_file;
  string top_file;
  string jobType;

  vector<AtomsRes> atoms;
 
  vector<string> itp_files;
  vector<string> isolabels;
  vector<int> chrom_Clist;
  vector<int> amideI_Clist;

  int namideI;
  int nchrom;
  int natoms;

  rvec *x;
  rvec box;

  void writeJ();
  void findAmideI();
  void amideIJob();

};

#endif
