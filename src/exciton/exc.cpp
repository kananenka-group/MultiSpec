#include "exc.h"

using namespace std;

Exc::Exc(string h_file_name, string d_file_name, int nchr, int nt, double deltaT,
         double corrT, double trel, int nv, double tsep, bool irs, bool ramans):
         Hfile(h_file_name), Dfile(d_file_name), nchrom(nchr), ntime(nt), dt(deltaT),
         tc(corrT), rlx_time(trel), navg(nv), sep_time(tsep), ir(irs), raman(ramans)

{
   // setting up some variables 
   int nfrmn;
   if(ir || raman){
      ncor = (int) (tc/dt);
      nsep = (int) (sep_time/dt); 
      nfrmn = navg*ncor+nsep*(navg-1);

      if(nfrmn > ntime){
        printf("\n Error! The input trajectory is too short, need %d more frames.\n\n",(nfrmn-ntime));
        exit(EXIT_FAILURE);
      }
   }

   printf("\n** Setting up calculation: **\n");
   printf("     Time step: %7.5f [ps] Correlation time: %7.5f [ps] \n", dt, tc);
   printf("     Statistical averaging will be performed using %d slices \n",navg);
   printf("     Slices will be separated by %7.5f [ps] \n",sep_time);
   printf("     Input number of frames: %d \n",ntime);
   printf("     %d frames are needed to calculate spectra \n",nfrmn);

   // Open files
   hinfile.open(Hfile, ios::binary);
   if(hinfile.fail()){
      printf(" Error! Could not open file: %s \n",Hfile.c_str());
      exit(EXIT_FAILURE);
   }

   if(ir) {
      dinfile.open(Dfile, ios::binary);
      if(dinfile.fail()){
         printf(" Error! Could not open file: %s \n",Dfile.c_str());
         exit(EXIT_FAILURE);
      }
   }

   if(raman) {
      pinfile.open(Pfile, ios::binary);
      if(pinfile.fail()){
         printf(" Error! Could not open file: %s \n",Pfile.c_str());
         exit(EXIT_FAILURE);
      }
   }

}

void Exc::run()
{

   // define some common variables
   ndim1   = nchrom*(nchrom+1)/2;
   nchrom2 = nchrom*nchrom;

   F.resize(nchrom2);
   H1.resize(nchrom2);

   // run jobs...
   if(ir) FTIR();
   if(raman) Raman();

}

void Exc::Raman()
{
   printf("\n** Raman module **\n");
   printf("     Relaxation time: %7.5f [ps] \n",rlx_time);

   plz_xx0.resize(nchrom);
   plz_xy0.resize(nchrom); 
   plz_xz0.resize(nchrom);
   plz_yy0.resize(nchrom);
   plz_yz0.resize(nchrom);
   plz_zz0.resize(nchrom);
   plz_xx.resize(nchrom);
   plz_xy.resize(nchrom); 
   plz_xz.resize(nchrom);
   plz_yy.resize(nchrom);
   plz_yz.resize(nchrom);
   plz_zz.resize(nchrom);

   VVT.resize(ncor);
   fill_n(VVT.begin(), ncor, complex_zero);
   VHT.resize(ncor);
   fill_n(VHT.begin(), ncor, complex_zero);

   VVw.resize(NFFT);
   fill_n(VVw.begin(), NFFT, complex_zero);
   VHw.resize(NFFT);
   fill_n(VHw.begin(), NFFT, complex_zero);

   printf("\n** Calculating VV and VH Raman TCFs **\n");
   int ist;
   for(int ii=0; ii<navg; ++ii){
      ist = ii*(ncor+nsep);
      fill_n(F.begin(), nchrom2, complex_zero);
      for(int jj=0; jj<nchrom;++jj)  F[jj*nchrom+jj] = complex_one;
      readHf(1);
      readPf(1,true);
      calcRm(0);
      for(int tt=1; tt<ncor; ++tt){
         readHf(1);
         readPf(1,false);
         moveF();
         calcRm(tt);
      }
      readHf(nsep);
      readPf(nsep,false);
   }

   fgrid1D();

   // FFT to get spectra
   FFT1D(VVT, VVw, dt, NFFT);
   FFT1D(VHT, VHw, dt, NFFT);

   // calculate depolarization ratio
   double ivv, ivh, dr;
   ivv = simpsonInt(wgrid1d, VVw);
   ivh = simpsonInt(wgrid1d, VHw);
   dr = ivh/ivv;
   printf("     Depolarization ratio %7.5f \n",dr);

   // print stuff..
   printTCF1D("vvt.dat", "VV", VVT);
   printTCF1D("vht.dat", "VH", VHT);
   printRamT();
   printRamS();
   printIw1D("vvw.dat", "VV", VVw);
   printIw1D("vhw.dat", "VH", VHw);

}

void Exc::FTIR()
{
//
   printf("\n** FTIR module **\n");
   printf("     Relaxation time: %7.5f [ps] \n",rlx_time);

   mu1_x.resize(nchrom);
   mu1_y.resize(nchrom);
   mu1_z.resize(nchrom);
   mu1_x0.resize(nchrom);
   mu1_y0.resize(nchrom);
   mu1_z0.resize(nchrom);

   mR1D.resize(ncor);
   fill_n(mR1D.begin(), ncor, complex_zero);

   IRw.resize(NFFT);
   fill_n(IRw.begin(), NFFT, complex_zero);

   printf("\n** Calculating dipole-dipole TCF **\n"); 
   int ist;
   for(int ii=0; ii<navg; ++ii){
      ist = ii*(ncor+nsep); 
      fill_n(F.begin(), nchrom2, complex_zero);
      for(int jj=0; jj<nchrom;++jj)  F[jj*nchrom+jj] = complex_one;
      readHf(1);
      readDf(1,true);
      calcR1D(0);
      for(int tt=1; tt<ncor; ++tt){
         readHf(1);
         readDf(1,false);
         moveF();
         calcR1D(tt);   
      }
      readHf(nsep);
      readDf(nsep,false);
   }

   fgrid1D();

   FFT1D(mR1D, IRw, dt, NFFT);

   printTCF1D("irt.dat", "IR", mR1D);
   printIw1D("irw.dat", "IR", IRw);


}

void Exc::calcR1D(int ti)
{
  MKL_INT lda;
  lda = (MKL_INT) nchrom;

  vector<complex<double>> mu0(nchrom, complex_zero);
  vector<complex<double>> mut(nchrom, complex_zero);
  vector<complex<double>> work(nchrom, complex_zero);

  complex<double> cx, cy, cz;

  // x
  for(int ii=0; ii<nchrom; ++ii) mu0[ii] = complex_one*mu1_x0[ii];
  for(int ii=0; ii<nchrom; ++ii) mut[ii] = conj(complex_one*mu1_x[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &mu0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &mut[0], 1, &work[0], 1, &cx); 

  // y
  for(int ii=0; ii<nchrom; ++ii) mu0[ii] = complex_one*mu1_y0[ii];
  for(int ii=0; ii<nchrom; ++ii) mut[ii] = conj(complex_one*mu1_y[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &mu0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &mut[0], 1, &work[0], 1, &cy); 

  // z
  for(int ii=0; ii<nchrom; ++ii) mu0[ii] = complex_one*mu1_z0[ii];
  for(int ii=0; ii<nchrom; ++ii) mut[ii] = conj(complex_one*mu1_z[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &mu0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &mut[0], 1, &work[0], 1, &cz); 

  double dtc = (double) ti;
  double exptc = exp(-1.0*dt*dtc/(2.0*rlx_time));
  complex<double> tmptcf = (cx+cy+cz)*exptc; ///3.0;
  mR1D[ti] += tmptcf;

}

void Exc::moveF()
{
  // 1. Diagonalize the Hamiltonian
  vector<double> W(nchrom,0.0);
  vector<double> Ht(nchrom2,0.0);

  MKL_INT info, lda;
  lda = (MKL_INT) nchrom;

  memcpy(&Ht[0], &H1[0], sizeof(double)*nchrom2);

  if (info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', lda, &Ht[0], lda, &W[0]) != 0)
  {
     printf("Error! LAPACKE_dsyeved returned \n");
     exit(EXIT_FAILURE);
  }

  // 2. Get exp(iE_nt/hbar) matrix from eigenvalues
  complex<double> arg;
  vector<complex<double>> U(nchrom2, complex_zero);

  for (int ii=0; ii<nchrom; ii++){
     arg = img*W[ii]*dt/HBAR;
     U[ii*nchrom+ii] = exp(arg);
   }

   // 3. Get expiH
   vector<complex<double>> evec(nchrom2, complex_zero);
   vector<complex<double>> work(nchrom2, complex_zero);
   vector<complex<double>> eiH(nchrom2, complex_zero); 

   for (int i=0; i<nchrom; i++ ){
      for (int j=0; j<nchrom; j++){
         evec[i*nchrom+j].real(Ht[i*nchrom+j]);
         evec[i*nchrom+j].imag(0.0);
      }
   }
   
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lda, lda, 
               lda, &complex_one, &evec[0], lda, &U[0], lda, &complex_zero, 
               &work[0], lda);

   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, lda, lda, 
               lda, &complex_one, &work[0], lda, &evec[0], lda, &complex_zero, 
               &eiH[0], lda);

   // 4. Propagate F
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lda, lda,
               lda, &complex_one, &eiH[0], lda, &F[0], lda, &complex_zero,
               &work[0], lda);

   memcpy(&F[0], &work[0], sizeof(complex<double>)*nchrom2);

}

void Exc::printTCF1D(string fname, string stype, vector<complex<double>> R1D)
{
   ofstream o_ir_file;
   o_ir_file.open(fname);
   if(!o_ir_file.is_open()){
      printf(" Error! Cannot open file: %s \n",fname.c_str());
      exit(EXIT_FAILURE); 
   }

   printf("     Writing %s TCF into %s \n",stype.c_str(),fname.c_str());
   for(int t1=0; t1<ncor; ++t1)
      o_ir_file << t1*dt << "  " << R1D[t1].real()/navg << "  " << R1D[t1].imag()/navg << endl;
    
   o_ir_file.close();
      
} 

void Exc::readHf(int nread)
{
   vector<float> Htmp;
   Htmp.resize(ndim1*nread);

   hinfile.read(reinterpret_cast<char*>(&Htmp[0]), Htmp.size()*sizeof(float));

   if(nread==1){
      int jj;
      for(int i=0; i<nchrom; ++i){
         jj = i*nchrom - i*(i+1)/2;
         for(int j=i; j<nchrom; j++){
            H1[i*nchrom + j] = (double) Htmp[jj+j];
            H1[j*nchrom + i] = H1[i*nchrom + j];
         }
      }
   }
 

}

void Exc::readDf(int nread, bool st)
{
   int size = 3*nchrom*nread;

   vector<float> Mtmp;
   Mtmp.resize(size);

   dinfile.read(reinterpret_cast<char*>(&Mtmp[0]), Mtmp.size()*sizeof(float));
 
   if(nread==1){
      for(int ii=0; ii<nchrom; ++ii){
         mu1_x[ii] = (double) Mtmp[3*ii];
         mu1_y[ii] = (double) Mtmp[3*ii+1];
         mu1_z[ii] = (double) Mtmp[3*ii+2];
      }

      if(st){
         memcpy(&mu1_x0[0], &mu1_x[0], sizeof(double)*nchrom);
         memcpy(&mu1_y0[0], &mu1_y[0], sizeof(double)*nchrom);
         memcpy(&mu1_z0[0], &mu1_z[0], sizeof(double)*nchrom);
      }
   }

}

void Exc::printIw1D(string specf, string stype, vector<complex<double>> Iw)
{
   ofstream o_iw_file;
   o_iw_file.open(specf);
   if(!o_iw_file.is_open()){
      printf(" Error! Cannot open file: %s \n",specf.c_str());
      exit(EXIT_FAILURE);
   }

   printf("     Writing %s Spectra into %s \n",stype.c_str(),specf.c_str());
   for(int i=NFFT/2, j=0; i<NFFT; ++i, j++)
      o_iw_file << wgrid1d[j] << "  " << Iw[i].real()/navg << endl;

   for(int i=0, j=NFFT/2; i<NFFT/2; ++i, ++j)
      o_iw_file << wgrid1d[j] << "  " << Iw[i].real()/navg << endl;

   o_iw_file.close();

}

void Exc::readPf(int nread, bool st)
{
   int size = 6*nchrom*nread;

   vector<float> Ptmp;
   Ptmp.resize(size);

   pinfile.read(reinterpret_cast<char*>(&Ptmp[0]), Ptmp.size()*sizeof(float));

   if(nread==1){
      for(int ii=0; ii<nchrom; ++ii){
         plz_xx[ii] = (double) Ptmp[6*ii];
         plz_xy[ii] = (double) Ptmp[6*ii+1];
         plz_xz[ii] = (double) Ptmp[6*ii+2];
         plz_yy[ii] = (double) Ptmp[6*ii+3];
         plz_yz[ii] = (double) Ptmp[6*ii+4];
         plz_zz[ii] = (double) Ptmp[6*ii+5];
      }

      if(st){
         memcpy(&plz_xx0[0], &plz_xx[0], sizeof(double)*nchrom);
         memcpy(&plz_xy0[0], &plz_xy[0], sizeof(double)*nchrom);
         memcpy(&plz_xz0[0], &plz_xz[0], sizeof(double)*nchrom);
         memcpy(&plz_yy0[0], &plz_yy[0], sizeof(double)*nchrom);
         memcpy(&plz_yz0[0], &plz_yz[0], sizeof(double)*nchrom);
         memcpy(&plz_zz0[0], &plz_zz[0], sizeof(double)*nchrom);
      }
   }

}

void Exc::calcRm(int ti)
{
  MKL_INT lda;
  lda = (MKL_INT) nchrom;

  vector<complex<double>> tp0(nchrom, complex_zero);
  vector<complex<double>> tpt(nchrom, complex_zero);
  vector<complex<double>> work(nchrom, complex_zero);

  complex<double> cxxxx, cyyyy, czzzz, iiii;
  complex<double> cxxyy, cyyxx, cyyzz, czzyy, cxxzz, czzxx, iijj;
  complex<double> cxyxy, cxzxz, cyzyz, ijij;

  // xxxx
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_xx0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_xx[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cxxxx); 

  // yyyy
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_yy0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_yy[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cyyyy); 

  // zzzz
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_zz0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_zz[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &czzzz); 

  iiii = cxxxx + cyyyy + czzzz;

  // xxyy
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_xx0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_yy[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cxxyy);

  // yyxx
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_yy0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_xx[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cyyxx);

  // yyzz
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_yy0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_zz[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cyyzz);

  // zzyy
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_zz0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_yy[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &czzyy);

  // xxzz
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_xx0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_zz[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cxxzz); 

  // zzxx
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_zz0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_xx[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &czzxx);

  iijj = cxxyy + cyyzz + cxxzz + cyyxx + czzyy + czzxx;

  // xyxy
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_xy0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_xy[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cxyxy);

  // xzxz
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_xz0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_xz[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cxzxz);

  // yzyz
  for(int ii=0; ii<nchrom; ++ii) tp0[ii] = complex_one*plz_yz0[ii];
  for(int ii=0; ii<nchrom; ++ii) tpt[ii] = conj(complex_one*plz_yz[ii]);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F[0], lda, &tp0[0], 1, &complex_zero, &work[0], 1);
  cblas_zdotu_sub(lda, &tpt[0], 1, &work[0], 1, &cyzyz);

  ijij = cxyxy + cxzxz + cyzyz;

  double dtc = (double) ti;
  double exptc = exp(-1.0*dt*dtc/(2.0*rlx_time));

  VVT[ti] += (3.0*iiii + iijj + 4.0*ijij)*exptc/15.0;
  VHT[ti] += (2.0*iiii - iijj + 6.0*ijij)*exptc/30.0;

}

void Exc::fgrid1D()
{
   wgrid1d.resize(NFFT);
   for(int i=NFFT/2, j=0; i<NFFT; ++i, ++j)
       wgrid1d[j] = 2*M_PI*HBAR*(i-NFFT)/(dt*NFFT);
   
   for(int i=0, j=NFFT/2; i<NFFT/2; ++i, ++j)
      wgrid1d[j] = 2*M_PI*HBAR*i/(dt*NFFT);

}

void Exc::printRamT()
{
   ofstream o_r_file;
   o_r_file.open("isot.dat");
   if(!o_r_file.is_open()){
      printf(" Error! Cannot open file: isot.dat \n");
      exit(EXIT_FAILURE);
   }

   printf("     Writing Raman isotropic TCF into isot.dat \n");
   for(int t1=0; t1<ncor; ++t1)
      o_r_file << t1*dt << "  " << (VVT[t1].real() - 0.75*VHT[t1].real())/navg << "  " << (VVT[t1].imag() - 0.75*VHT[t1].imag())/navg << endl;

   o_r_file.close();

   ofstream o_ru_file;
   o_ru_file.open("unpt.dat");
   if(!o_ru_file.is_open()){
      printf(" Error! Cannot open file: unpt.dat \n");
      exit(EXIT_FAILURE);
   }

   printf("     Writing Raman unpolarized TCF into unpt.dat \n");
   for(int t1=0; t1<ncor; ++t1)
      o_ru_file << t1*dt << "  " << (VVT[t1].real() + VHT[t1].real())/navg << "  " << (VVT[t1].imag() + VHT[t1].imag())/navg << endl;

   o_ru_file.close();

}

void Exc::printRamS()
{
   ofstream o_isw_file;
   o_isw_file.open("isow.dat");
   if(!o_isw_file.is_open()){
      printf(" Error! Cannot open file: isow.dat \n");
      exit(EXIT_FAILURE);
   }

   printf("     Writing Raman isotropic Spectra into isow.dat \n");
   for(int i=NFFT/2, j=0; i<NFFT; ++i, j++)
      o_isw_file << wgrid1d[j] << "  " << (VVw[i].real() - (4.0/3.0)*VHw[i].real())/navg << endl;

   for(int i=0, j=NFFT/2; i<NFFT/2; ++i, ++j)
      o_isw_file << wgrid1d[j] << "  " << (VVw[i].real() - (4.0/3.0)*VHw[i].real())/navg << endl;

   o_isw_file.close();

   ofstream o_upw_file;
   o_upw_file.open("unpw.dat");
   if(!o_upw_file.is_open()){
      printf(" Error! Cannot open file: unpw.dat \n");
      exit(EXIT_FAILURE);
   }

   printf("     Writing Raman unpolarized Spectra into unpw.dat \n");
   for(int i=NFFT/2, j=0; i<NFFT; ++i, j++)
      o_upw_file << wgrid1d[j] << "  " << (VVw[i].real() + VHw[i].real())/navg << endl;

   for(int i=0, j=NFFT/2; i<NFFT/2; ++i, ++j)
      o_upw_file << wgrid1d[j] << "  " << (VVw[i].real() + VHw[i].real())/navg << endl;

   o_upw_file.close();

   
}
