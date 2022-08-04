#ifndef WM_H
#define WM_H

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm> // for transform
#include <iterator> // begin, end, ostream_iterator
#include <numeric>  // iota
#include "wmaps.h"
#include "wmapb.h"
#include "../traj/traj.h"
#include "../vec/vecManip.h"
#include "../const/const.h"
#include "../util/util.h"

using namespace std;

class water{

public:
   water(string, string, string, string, string, string, 
         string, int, int, bool, bool, bool, bool, bool,
         bool, float, float);

  ~water();

private:
  string water_model_name, water_model_name_inp;
  string jobType;
  string traj_file;
  string gro_file;
  string ams_file;

  WmapS wms;
  WmapB wmb;

  int nframes;
  int nd2o = 0;
  
  bool ir;
  bool raman;
  bool sfg;
  bool DoDv = false;
  bool intrac;
  bool intermcs;

  float tdSFG;
  float fc;

  int atoms_in_mol;
  int natoms;
  int nwater;
  int nchroms = 0;
  int nchromb = 0;
  int nchromt;
  int nchromt2;
  int ndim;
  int offs = 0;

  rvec *x;
  rvec box;
  rvec *vOHa, *vOHu, *eft; 

  vector<int> oxyInd;
  vector<int> mH2O;
  vector<int> mD2O;
  vector<int> mHOD;
  vector<int> woxyT;

  vector<float> ROH;

  vector<float> uChg;
  vector<float> aChg;
  vector<float> uMas;
  vector<float> aMas;

  vector<float> efs;
  vector<float> efb;
  vector<float> hf;
  vector<float> tdmuf;
  vector<float> plzbf;
  vector<float> fzf;

  vector<float> w10;
  vector<float> x10;
  vector<float> p10;
  vector<float> m10;
  vector<float> w20b;

  vector<float> omegas;
  vector<int> w_dist;

  int now, nhw, nmw;
  int nh2o, nhod;
  int nbins = 100;

  bool ws = true;
  bool wb = false;
  bool wf = false;;
  bool uncs = false;
  bool moveM=false;

  float total_charge;
  float total_mass;
  float inv_total_mass;
  float trdip;
  float pz;
  float rom;

  string water_map_name;
  string chromType;
  string w_dist_fname = "freq_hist.dat";

  ofstream houtfile;
  ofstream doutfile;
  ofstream poutfile;
  ofstream jobfile;
  ofstream fzoutfile;
  ofstream isofile;
  ofstream inpfile;

  vector<string> uAtoms;
  vector<string> aAtoms;

  void waterModel();
  void waterJob();
  void readGro();
  void readCharges();
  void IsoMix();
  void calcEf();
  void updateEx();
  void calcWXPM();
  double waterTDC(const rvec&, const rvec&, const rvec&, const rvec&, const rvec&);
  void com(rvec&);
  void intermC();
  void trDip();
  void trPol();
  void CalcSQuant();
  void moveMsite();
  void freqDist();

  void writeH();
  void writeD();
  void writeP();
  void writeFz();

};

#endif
