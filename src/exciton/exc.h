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
#include "ir2d.h"
#include "../const/const.h"
#include "../fft/dofft.h"
#include "../math/nummath.h"
#include "../util/util.h"

using namespace std;

class Exc{

public:
   Exc(string, string, string, int, int, int, int,
       double, double, double, double, double,
       double, double, double, double, bool, 
       bool, bool, bool, bool);
  ~Exc(){};

   void run();

private:
   void readHf(int);
   void readHs(int);
   void readDs(int);
   void readDf(int,bool);
   void readPf(int,bool);
   void readFzf(int,bool);
   void FTIR();
   void calc2DIR();
   void start2DIR();
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
   void buildM21();
   void expH(vector<complex<double>>&, vector<double>, int);
   void prop(vector<complex<double>>&, vector<double>, int, int);
   void assignD(int);
   complex<double> calcR2D_R1(int);
   complex<double> calcR2D_R2(int);
   void save_F1_t0t1_mu0(int);
   void assignDpol(vector<complex<double>> &, vector<complex<double>> &,
                   vector<complex<double>> &, vector<complex<double>> &,
                   vector<complex<double>> &, vector<complex<double>> &,
                   int, int);
   void writeR2D();
   void write2DRabs();
   void write2Dout(vector<complex<double>>, string, string, int);

   vector<double> H1;
   vector<double> mu1_x;
   vector<double> mu1_y;
   vector<double> mu1_z;
   vector<double> mu01_x;
   vector<double> mu01_y;
   vector<double> mu01_z;
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

   // 2D 
   vector<double> H2;                   // two-exciton Hamiltonian read from Hfile
   vector<complex<double>> F1_t0t1;     // one-exciton Hamiltonian propagated t0 -> t1
   vector<complex<double>> F1_t1t2;     // one-exciton Hamiltonian propagated t1 -> t2
   vector<complex<double>> F1_t2t3;     // one-exciton Hamiltonian propagated t2 -> t3 
   vector<complex<double>> F2_t1t2;     // two-exciton Hamiltonian propagated t1 -> t2
   vector<complex<double>> F2_t2t3;     // two-exciton Hamiltonian propagated t2 -> t3
   vector<complex<double>> R2D_R1;
   vector<complex<double>> R2D_R2;
   vector<complex<double>> R1D;
   vector<complex<double>> F1_t0t1_mu01_0_x;
   vector<complex<double>> F1_t0t1_mu01_0_y;
   vector<complex<double>> F1_t0t1_mu01_0_z;
   vector<double> mu12_x;
   vector<double> mu12_y;
   vector<double> mu12_z;
   vector<double> mu01_t0_x;
   vector<double> mu01_t0_y;
   vector<double> mu01_t0_z;
   vector<double> mu01_t1_x;
   vector<double> mu01_t1_y;
   vector<double> mu01_t1_z;
   vector<double> mu01_t2_x;
   vector<double> mu01_t2_y;
   vector<double> mu01_t2_z;
   vector<double> mu01_t3_x;
   vector<double> mu01_t3_y;
   vector<double> mu01_t3_z;
   vector<double> mu12_t2_x;
   vector<double> mu12_t2_y;
   vector<double> mu12_t2_z;
   vector<double> mu12_t3_x;
   vector<double> mu12_t3_y;
   vector<double> mu12_t3_z;
   vector<complex<double>> mu0_eg;
   vector<complex<double>> mu1_eg;
   vector<complex<double>> mu2_eg;
   vector<complex<double>> mu3_eg;
   vector<complex<double>> mu2_ce;
   vector<complex<double>> mu3_ce;

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
   int n1ex;
   int n1ex2;
   int n2ex;
   int n2ex2;
   int ntime;
   int nchrom2;
   int ncor;
   int navg;
   int nstart;
   int nsep;
   int ndim1;
   int NFFT;
   int nt1t3;
   int nt2;
   int nfrmn=0;

   double dt;
   double tc;
   double T1_rlx;
   double sep_time;
   double start_time;
   double w_avg;
   double anharm;
   double t1t3;
   double t2;

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
