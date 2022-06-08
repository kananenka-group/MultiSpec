#ifndef INPUT_H
#define INPUT_H

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

  double getdt()   const { return dt; }
  double gettc()   const { return tc; }
  double getTrlx() const { return trlx; }
  double getTsep() const { return tsep; }

  bool ifIR()     const { return ir; }
  bool ifRaman()  const { return raman; }
  bool ifSFG()    const { return sfg; }

  string getHfile() const { return Hfile; }
  string getDfile() const { return Dfile; }
  string getPfile() const { return Pfile; }

private:
  int nframe = 1;
  int ntime = 1;
  int nchrom = 1;
  int navg = 1;
 
  double dt = 0.01;  // time step between frames in ps
  double tc = 1.00;  // correlation time in ps
  double trlx = 0.2;
  double tsep = 2.0; // separation time in ps

  bool ir = false;     // calc IR spectra
  bool raman = false;  // calc Raman spectra
  bool sfg = false;    // calc SFG spectra

  string Hfile = "Hamiltonian.bin";
  string Dfile = "Dipole.bin";
  string Pfile = "Polarizability.bin";


};

#endif
