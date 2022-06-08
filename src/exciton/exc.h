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

using namespace std;

class Exc{

public:
   Exc(string, string, int, int, double, double, double, int, double, 
       bool, bool);
  ~Exc(){};

   void run();

private:
   void readHf(int);
   void readDf(int,bool);
   void readPf(int,bool);
   void FTIR();
   void Raman();
   void moveF();
   void calcR1D(int);
   void calcRm(int);
   void printTCF1D(string, string, vector<complex<double>>);
   void printIw1D(string, string, vector<complex<double>>);
   void fgrid1D();
   void printRamT();
   void printRamS();

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

   vector<complex<double>> F;
   vector<complex<double>> mR1D;
   vector<complex<double>> VVT;
   vector<complex<double>> VHT;
   vector<complex<double>> IRw;
   vector<complex<double>> VVw;
   vector<complex<double>> VHw;

   double tc;
   double dt;
   double rlx_time;
   double sep_time;

   int nchrom;
   int ntime;
   int nchrom2;
   int ncor;
   int navg;
   int nsep;
   int ndim1;
   int NFFT = 4096;

   bool ir;
   bool raman;

   // complex constants
   const complex<double> img          = {0.,1.};
   const complex<double> complex_one  = {1.,0.};
   const complex<double> complex_zero = {0.,0.};

   string Hfile = "Hamiltonian.bin";
   string Dfile = "Dipole.bin";
   string Pfile = "Polarizability.bin";

   ifstream hinfile;
   ifstream dinfile;
   ifstream pinfile;

};

#endif
