#ifndef WM_H
#define WM_H

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
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
   water(string, int, string, string, string, string, string, 
         string, bool, bool, bool, int, float, bool, float);

  ~water();

private:
  WmapS wms;
  WmapB wmb;

  int atoms_in_mol;
  int natoms;
  int nframes;
  int nwater;
  int nchroms = 0;
  int nchromb = 0;
  int nchromt;
  int nchromt2;
  int ndim;
  int offs = 0;

  const rvec *x;
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

  int now, nhw, nmw;
  int nh2o, nhod;
  int nd2o = 0;

  bool ws = true;
  bool wb = false;
  bool wf = false;;
  bool uncs = false;
  bool DoDv = false;

  bool ir;
  bool raman;
  bool sfg;

  float total_charge;
  float total_mass;
  float inv_total_mass;
  float trdip;
  float pz;
  float fc;
  float tdSFG;

  string water_model_name, water_model_name_caps;
  string water_map_name;
  string traj_file;
  string chromType;
  string jobType;
  string gro_file;
  string ams_file;

  ofstream houtfile;
  ofstream doutfile;
  ofstream poutfile;
  ofstream jobfile;
  ofstream fzoutfile;
  ofstream isofile;

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

  void writeH();
  void writeD();
  void writeP();
  void writeFz();

};

#endif
