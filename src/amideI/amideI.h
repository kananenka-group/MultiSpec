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
   amideI(string, string, vector<string>);

  ~amideI();

private:

  System s;

  string gro_file;
  string traj_file;
 
  vector<string> itp_files;

  int namideI;

  rvec *x;
  rvec box;

  vector<int> aAmideI;
  vector<int> residue_numbers;
  vector<int> res_start_numb;
  vector<string> atoms;
  vector<string> residue_names;
  vector<vector<string>> resAmideI;
  vector<vector<int>> resNAmideI;

  void writeJ();

};

#endif
