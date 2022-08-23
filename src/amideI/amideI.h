#ifndef AMIDE_H
#define AMIDE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <algorithm> // for transform
#include <iterator> // begin, end, ostream_iterator
#include <numeric>  // iota
#include "nnmap.h"
#include "../system/system.h"
#include "../traj/traj.h"
#include "../vec/vecManip.h"
#include "../const/const.h"
#include "../util/util.h"

using namespace std;

class amideI{

public:
   amideI(string, string, vector<string>, string, string,
          vector<string>, string, int, int, bool);

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

  NNmap map_nn;

  int nframes;
  int startframe=1;

  bool ir;

  ofstream houtfile;
  ofstream doutfile;
  ofstream jobfile;

  vector<float> hf;
  vector<float> tdmuf;
  vector<float> diag_w_nn;
  vector<float> diag_w_e;
  //vector<float> phi;
  //vector<float> psi;

  //vector<int> angleID;

  int namideI;
  int nchrom;
  int nchrom2;
  int natoms;

  bool amdsiso = false;

  rvec *x;
  rvec box;

  void writeJ();
  void findAmideI();
  void amideIJob();
  void nnfs();
  void CalcSQuant();
  void calcAngles(const int , const int , const rvec *, float &, float &);
  float calcDihedral(const rvec &, const rvec &,
                     const rvec &, const rvec &);
  uint search2(const int , const int , const string &, const int);
  float calc_N_NN_shift(const int , const int , const rvec *);
  float calc_C_NN_shift(const int , const int , const rvec *);

};

#endif
