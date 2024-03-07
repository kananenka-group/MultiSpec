#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class Input{

public:
   Input(int argc, char ** argv);
  ~Input() {};

  int getNchrom() const { return nchrom; }
  int getNtime()  const { return ntime;  }
  int getNfrm()   const { return nframe; }
  int getNavg()   const { return navg; }
  int getNFFT()   const { return NFFT; }

  double getdt()   const { return dt; }
  double gettc()   const { return tc; }
  double getts()   const { return ts; }
  double getTrlx() const { return trlx; }
  double getTsep() const { return tsep; }
  double getWavg() const { return w_avg; }
  double getAnh()  const { return anharm; }
  double gett1t3() const { return t1t3; }
  double gett2()   const { return t2; }

  bool ifIR()     const { return ir; }
  bool ifIR2D()   const { return ir2d; }
  bool ifRaman()  const { return raman; }
  bool ifSFG()    const { return sfg; }
  bool ifSD()     const { return sd; }

  string getHfile() const { return Hfile; }
  string getDfile() const { return Dfile; }
  string getPfile() const { return Pfile; }

private:
  int nframe = 1;
  int ntime = 1;
  int nchrom = 1;
  int navg = 1;
  int NFFT = 4096;
 
  double dt = 0.01;  // time step between frames in ps
  double tc = 1.00;  // correlation time in ps
  double trlx = 0.2;
  double tsep = 2.0; // separation time in ps
  double w_avg = 0.0;
  double ts = 0.0;
  double anharm = 0.0; // diagonal anharmonicity for 2d ir
  double t1t3 = 2.5;   // t1 and t3 time in 2D IR in ps
  double t2   = 0.0;   // t2 time in 2D IR

  bool ir = false;     // calc IR spectra
  bool ir2d = false;  // calculate 2D IR spectra
  bool raman = false;  // calc Raman spectra
  bool sfg = false;    // calc SFG spectra
  bool sd = false;

  string Hfile = "Hamiltonian.bin";
  string Dfile = "Dipole.bin";
  string Pfile = "Polarizability.bin";
  string Jfile = "job.bin";

  ifstream jobfile;
  ofstream inpfile;


};

#endif
