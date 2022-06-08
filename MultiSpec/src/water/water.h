#ifndef WM_H
#define WM_H

#include <iostream>
#include <sstream>
#include <string>
#include "wmaps.h"
#include "../traj/traj.h"
#include "../vec/vecManip.h"
#include "../const/const.h"

using namespace std;

class water{

public:
   water(string, string, int, string, string, string, string, string,
         bool, bool, bool);
  ~water() {};

  int getNwatoms()       const { return atoms_in_mol; }
  string getWModelName() const { return water_model_name_caps; }

private:
  Wmap wm;

  int atoms_in_mol, atoms_in_mol_wm;
  int natoms;
  int nframes;
  int nwater;
  int nchrom;
  int ndim;

  vector<int> oxyInd;
  vector<float> uChg;
  vector<float> aChg;
  vector<float> ef;
  vector<float> hf;
  vector<float> ht;
  vector<float> tdmut;
  vector<float> tdmuf;
  vector<float> plzbt;
  vector<float> plzbf;

  int now, nhw, nmw;

  bool ws;
  bool wb;
  bool wf;
  bool pure;
  bool ir;
  bool raman;
  bool sfg;

  float total_charge;
  float trdip;
  float pz;

  string water_model_name, water_model_name_caps;
  string OH_water_map_name;
  string traj_file;
  string chromType;
  string jobType;
  string gro_file;
  string chg_file;

  vector<string> uAtoms;
  vector<string> aAtoms;

  void waterModel();
  void waterChrom();
  void waterJob();
  void readGro();
  void readCharges();
  void updateEx();

};

#endif
