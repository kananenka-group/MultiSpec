#ifndef EXC_H
#define EXC_H

#include <algorithm>
#include <iterator>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include "mkl.h"
#include <omp.h>
#include <pthread.h>
#include "../const/const.h"
#include "../fft/dofft.h"
#include "../math/nummath.h"

using namespace std;

class Exc{

public:
   Exc(string, string, string, int, int, int, 
       double, double, double, double, double,
       double, double, bool, bool, bool, bool, 
       bool);
  ~Exc(){};

   void run();

private:
   void readHf(int);
   void readDf(int,bool);
   void readPf(int,bool);
   void readFzf(int,bool);
   void FTIR();
   void calc2DIR();
   void Raman();
   void SFG();
   void moveF();
   void calcR1D(int);
   void calcRm(int);
   void calcSFG(int);
   void printTCF1D(string, string, vector<complex<double>>);
   void printIw1D(string, string, vector<complex<double>>, int);
   void printSd1D(string, string, vector<double>);
   void fgrid1D();
   void printRamT();
   void printRamS();
   void scaleTDM();
   void buildH2();
   void sdIR();
   void sdSFG();
   void sdRaman();
   int get2nx(int, int);

   vector<double> H1;
   vector<double> mu1_x;
   vector<double> mu1_y;
   vector<double> mu1_z;
   vector<double> mu1_x0;
   vector<double> mu1_y0;
   vector<double> mu1_z0;
   vector<double> plz_xx;
   vector<double> plz_xy;
   vector<double> plz_xz;
   vector<double> plz_yy;
   vector<double> plz_yz;
   vector<double> plz_zz;
   vector<double> plz_xx0;
   vector<double> plz_xy0;
   vector<double> plz_xz0;
   vector<double> plz_yy0;
   vector<double> plz_yz0;
   vector<double> plz_zz0;
   vector<double> wgrid1d;
   vector<double> fzt;
   vector<double> fzt0;
   vector<double> evecsr;
   vector<double> evals;
   vector<double> sdr_ir;
   vector<double> sdr_sfg;
   vector<double> sdr_vv;
   vector<double> sdr_vh;
   vector<double> sdr_w;

   vector<double> H2;

   vector<complex<double>> F;
   vector<complex<double>> mR1D;
   vector<complex<double>> VVT;
   vector<complex<double>> VHT;
   vector<complex<double>> IRw;
   vector<complex<double>> VVw;
   vector<complex<double>> VHw;
   vector<complex<double>> sspt;
   vector<complex<double>> pppt;
   vector<complex<double>> sspw;
   vector<complex<double>> pppw;
   vector<complex<double>> yyzt;
   vector<complex<double>> yyzw;
   vector<complex<double>> spst;
   vector<complex<double>> spsw;

   string Hfile  = "Hamiltonian.bin";
   string Dfile  = "Dipole.bin";
   string Pfile  = "Polarizability.bin";
   string Fzfile = "Fz.bin";

   int nchrom;
   int n2ex;
   int n2ex2;
   int ntime;
   int nchrom2;
   int ncor;
   int navg;
   int nstart;
   int nsep;
   int ndim1;
   int NFFT = 4096;

   double dt;
   double tc;
   double T1;
   double T2;
   double sep_time;
   double start_time;
   double w_avg;
   double anharm;

   bool ir;
   bool ir2d;
   bool raman;
   bool sfg;
   bool sd;

   // complex constants
   const complex<double> img          = {0.,1.};
   const complex<double> complex_one  = {1.,0.};
   const complex<double> complex_zero = {0.,0.};

   // real constants
   const double done  = 1.0;
   const double dzero = 0.0;

   ifstream hinfile;
   ifstream dinfile;
   ifstream pinfile;
   ifstream finfile;

};

#endif
