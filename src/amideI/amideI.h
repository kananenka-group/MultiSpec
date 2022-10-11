#ifndef AMIDE_H
#define AMIDE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <algorithm> 
#include <iterator> 
#include <numeric>  
#include "nnmap.h"
#include "nncmap.h"
#include "elmap.h"
#include "../system/system.h"
#include "../traj/traj.h"
#include "../vec/vecManip.h"
#include "../const/const.h"
#include "../util/util.h"

using namespace std;

class amideI{

public:
   amideI(string, string, string, vector<string>, string, string,
          vector<string>, string, string, string, int, int, 
          bool, float);

  ~amideI();

private:

  System s;
  
  string gro_file;
  string traj_file;
  string mapsFile;
  string top_file;
  string jobType;

  vector<AtomsRes> atoms;
 
  vector<string> itp_files;
  vector<string> isolabels;
  vector<int> chrom_Clist;
  vector<int> amideI_Clist;

  NNmap map_nn;
  NNCmap map_nnc;
  ELmap map_el;

  int nframes;
  int startframe=1;

  bool ir;

  float isoShift;

  ofstream houtfile;
  ofstream doutfile;
  ofstream jobfile;

  vector<float> hf;
  vector<float> tdmuf;
  vector<float> diag_w_nn;
  vector<float> diag_w_el;

  vector<int> chgSt;
  vector<int> grpInd;
  vector<float> grpQ;
  vector<vector<int> > exclude_list, exclude_group_list;

  rvec *vecCO, *vecCN, *trdipL, *dip;
  rvec *cog;
 
  //vector<float> phi;
  //vector<float> psi;

  //vector<int> angleID;

  int namideI;
  int nchrom;
  int nchrom2;
  int natoms;
  int ncog;

  bool amdsiso = false;
  bool save = false;

  rvec *x;
  rvec box;

  void writeJ();
  void writeH();
  void writeD();
  void findAmideI();
  void amideIJob();
  void nnfs();
  void elst();
  void elstCG();
  void couplings();
  void trDip();
  void updateExd();
  void CalcSQuant();
  vector<int> get_excludes(int, int);
  void calcAngles(const int , const int , const rvec *, float &, float &);
  float calcDihedral(const rvec &, const rvec &,
                     const rvec &, const rvec &);
  uint search2(const int , const int , const string &, const int);
  float calc_N_NN_CS(const int , const int , const rvec *, int);
  float calc_C_NN_CS(const int , const int , const rvec *, int);
  void getCCG(); 

};

#endif
