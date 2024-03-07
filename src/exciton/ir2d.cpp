#include "ir2d.h"

void Exc::start2DIR()
{
// perform some necessary checks...
  nstart = (int) (start_time/dt);
  nsep   = (int) (sep_time/dt);
  nt2    = (int) (t2/dt);
  nt1t3  = (int) (t1t3/dt); 

  // total number of frames needed
  nfrmn = nstart + (navg-1)*nsep + navg*(2*nt1t3+nt2);

  // check if we have enough frames
  if(nfrmn > ntime){
     printf("\n Error! Input trajectory is too short, need %d more frames.\n\n",(nfrmn-ntime));
     exit(EXIT_FAILURE);
  }

  // all is good print some info
  printf("\n** Setting up calculation: **\n");
  printf("   Time step: %7.5f [ps] \n",dt);
  printf("   t1 and t3 time: %7.5f [ps] \n",t1t3);
  printf("   t2 time: %7.5f [ps] \n",t2);
  printf("   Statistical averaging will be performed using %d slices \n",navg);
  printf("   Starting time : %7.5f [ps] \n",start_time);
  printf("   Slices will be separated by %7.5f [ps] \n",sep_time);
  printf("   %d frames are needed to calculate spectra \n",nfrmn);
  printf("   Input number of frames: %d \n",ntime);
}

void Exc::calc2DIR()
{
   printf("\n** 2D IR module **\n");

   // set up some variables
   // one-exciton subspace
   n1ex = nchrom;
   n1ex2= n1ex*n1ex;

   // two-exciton subspace
   n2ex = nchrom*(nchrom+1)/2;
   n2ex2= n2ex*n2ex;

   // check Hamiltonian input file
   hinfile.open(Hfile, ios::binary);
   if(hinfile.fail()){
       printf(" Error! Could not open file: %s \n",Hfile.c_str());
       exit(EXIT_FAILURE);
   }
   // don't need next 2 lines
   hinfile.clear();
   hinfile.seekg(0, ios::beg);
  
   // check Dipole input file 
   dinfile.open(Dfile, ios::binary);
   if(dinfile.fail()){
      printf(" Error! Could not open file: %s \n",Dfile.c_str());
      exit(EXIT_FAILURE);
   }
   // don't need next 2 lines
   dinfile.clear();
   dinfile.seekg(0, ios::beg);
   
   // relaxation time:
   printf("   Relaxation time T1: %5.3f [ps] \n",T1_rlx);
   printf("   Diagonal anharmonicity: %7.2f [cm-1]\n",anharm);
   
   // main loop...
   printf("\n** Calculating 2D IR response functions **\n");

   // resize some variables here...
   mu01_x.resize(n1ex, 0.0);
   mu01_y.resize(n1ex, 0.0);
   mu01_z.resize(n1ex, 0.0);
   mu12_x.resize(n1ex*n2ex, 0.0);
   mu12_y.resize(n1ex*n2ex, 0.0);
   mu12_z.resize(n1ex*n2ex, 0.0);
   mu01_t0_x.resize(n1ex, 0.0);
   mu01_t0_y.resize(n1ex, 0.0);
   mu01_t0_z.resize(n1ex, 0.0);
   mu01_t1_x.resize(n1ex, 0.0);
   mu01_t1_y.resize(n1ex, 0.0);
   mu01_t1_z.resize(n1ex, 0.0);
   mu01_t2_x.resize(n1ex, 0.0);
   mu01_t2_y.resize(n1ex, 0.0);
   mu01_t2_z.resize(n1ex, 0.0);
   mu01_t3_x.resize(n1ex, 0.0);
   mu01_t3_y.resize(n1ex, 0.0);
   mu01_t3_z.resize(n1ex, 0.0);
   mu12_t2_x.resize(n1ex*n2ex, 0.0);
   mu12_t2_y.resize(n1ex*n2ex, 0.0);
   mu12_t2_z.resize(n1ex*n2ex, 0.0);
   mu12_t3_x.resize(n1ex*n2ex, 0.0);
   mu12_t3_y.resize(n1ex*n2ex, 0.0);
   mu12_t3_z.resize(n1ex*n2ex, 0.0);

   F1_t0t1.resize(n1ex2, complex_zero);
   F1_t1t2.resize(n1ex2, complex_zero);
   F1_t2t3.resize(n1ex2, complex_zero);
   F2_t1t2.resize(n2ex2, complex_zero);
   F2_t2t3.resize(n2ex2, complex_zero);
   H1.resize(n1ex2, 0.0);
   H2.resize(n2ex2, 0.0);

   F1_t0t1_mu01_0_x.resize(n1ex*nt1t3, complex_zero);
   F1_t0t1_mu01_0_y.resize(n1ex*nt1t3, complex_zero);
   F1_t0t1_mu01_0_z.resize(n1ex*nt1t3, complex_zero);

   mu0_eg.resize(n1ex, complex_zero);
   mu1_eg.resize(n1ex, complex_zero);
   mu2_eg.resize(n1ex, complex_zero);
   mu3_eg.resize(n1ex, complex_zero);
   mu2_ce.resize(n1ex*n2ex, complex_zero);
   mu3_ce.resize(n1ex*n2ex, complex_zero);

   R1D.resize(nt1t3, complex_zero);
   R2D_R1.resize(nt1t3*nt1t3, complex_zero);
   R2D_R2.resize(nt1t3*nt1t3, complex_zero);

   int it0_min, it1_max, it2_max, it3_max;
   int it0, it1, it3;
   int ndx;

   for(int ii=0; ii<navg; ++ii){
      it0_min = nstart + ii*nsep;
      it1_max = it0_min + nt1t3-1;
      it2_max = it1_max + nt2;
      it3_max = it2_max + nt1t3 - 1;
      
      // begining new sample, reset F1
      fill_n(F1_t0t1.begin(), n1ex2, complex_zero);
      for(int i=0; i<n1ex; i++) F1_t0t1[i*n1ex+i] = complex_one;
      
      // propagate F1_t0t1 backward:
      for(it0 = it1_max; it0 >= it0_min; it0--){
         readDs(it0);
         assignD(0);
         save_F1_t0t1_mu0(it1_max-it0);
         if(it0 == it0_min) continue;
         readHs(it0);
         prop(F1_t0t1, H1, n1ex, -1);
      }

      // t1
      readDs(it1_max);
      assignD(1);
      //linear response

      // t1 -> t2
      fill_n(F1_t1t2.begin(), n1ex2, complex_zero);
      for(int i=0; i<n1ex; i++) F1_t1t2[i*n1ex+i] = complex_one;
      fill_n(F2_t1t2.begin(), n2ex2, complex_zero);
      for(int i=0; i<n2ex; i++) F2_t1t2[i*n2ex+i] = complex_one;

      for(int it2 = it1_max; it2<it2_max; it2++){
         readHs(it2);
         prop(F1_t1t2, H1, n1ex, 1);
         prop(F2_t1t2, H2, n2ex, 1);
      }
     
      // t2
      readDs(it2_max);
      assignD(2);

      // loop over t3
      fill_n(F1_t2t3.begin(), n1ex2, complex_zero);
      for(int i=0; i<n1ex; i++) F1_t2t3[i*n1ex+i] = complex_one;
      fill_n(F2_t2t3.begin(), n2ex2, complex_zero);
      for(int i=0; i<n2ex; i++) F2_t2t3[i*n2ex+i] = complex_one;
    
      for(it3 = it2_max; it3 <=it3_max; it3++){
         readDs(it3);
         assignD(3);
         // RDF here
         //#pragma omp parallel for private(ndx)
         for(it0=0; it0<nt1t3; it0++){
            ndx = it0*nt1t3 + it3 - it2_max;
            R2D_R1[ndx] += calcR2D_R1(it0);
            R2D_R2[ndx] += calcR2D_R2(it0); 
         }
         if(it3 == it3_max) continue; 
         readHs(it3);
         prop(F1_t2t3, H1, n1ex, 1);
         prop(F2_t2t3, H2, n2ex, 1);
      }

   }

   // account for dephasing
   printf(" !!!! Check it2 in dephasing, not clear \n");
   for(it1 = 0; it1 < nt1t3; ++it1){
      // linear
      R1D[it1] *= exp(-it1*dt/(2.0*T1_rlx))/(1.0*navg);
      for(it3 = 0; it3 < nt1t3; ++it3){
         ndx = it1*nt1t3 + it3;
         // rephasing
         R2D_R1[ndx] *= exp(-(it1+2.0*nt2+it3)*dt/(2.0*T1_rlx))/(1.0*navg);
         // non-rephasing
         R2D_R2[ndx] *= exp(-(it1+2.0*nt2+it3)*dt/(2.0*T1_rlx))/(1.0*navg);
      }
   }

   // finish and print spectra
   writeR2D();

   // FFT
   printf(" running FFTs\n");
   write2DRabs();

}

void Exc::expH(vector<complex<double>> &eH, vector<double> H, int N)
{
//
// Build exp(-iHt/hbar)
//
// 1. Diagonalize the Hamiltonian
   int N2 = N*N;
   int i, j;

   complex<double> arg; 
 
   vector<double> W(N,0.0);
   vector<double> Ht(N2,0.0);
   vector<complex<double>> evec(N2, complex_zero);
   vector<complex<double>> work(N2, complex_zero);

   MKL_INT info, lda;
   lda = (MKL_INT) N;

   memcpy(&Ht[0], &H[0], sizeof(double)*N2);
   
   info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', lda, &Ht[0], lda, &W[0]);
   if(info !=0){
      printf("Error! LAPACKE_dsyeved returned %d in Exc::expH \n",info);
      exit(EXIT_FAILURE);
   }

   fill_n(eH.begin(), N2, complex_zero);
  
   for(i=0; i<N; ++i){
      arg = -img*W[i]*dt/constants::HBAR;
      eH[i*N+i] = exp(arg);
   }

   for(i=0; i<N; i++){
      for(j=0; j<N; j++){
         evec[i*N+j].real(Ht[i*N+j]);
         evec[i*N+j].imag(0.0);
      }
   }

   // convert back to original basis
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lda, lda,
               lda, &complex_one, &evec[0], lda, &eH[0], lda, &complex_zero,
               &work[0], lda);

   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, lda, lda,
               lda, &complex_one, &work[0], lda, &evec[0], lda, &complex_zero,
               &eH[0], lda);
   
}

void Exc::prop(vector<complex<double>> &F, vector<double> H, int N, int dir)
{
//
// Propagate F matrix here
// dir = 1 -> forward propagation:   F(t+dt) = exp(-iHdt/hbar)*F(t)
// dir = -1 -> backward propagation: F(t-dt) = F(t)*exp(-iHdt/hbar)
//
    int N2 = N*N;

    vector<complex<double>> eiH(N2, complex_zero);
    vector<complex<double>> W(N2, complex_zero);
   
    expH(eiH,H,N);

    MKL_INT lda;
    lda = (MKL_INT) N;

    if(dir==1){
       // forward propagator
       cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lda, lda,
                   lda, &complex_one, &eiH[0], lda, &F[0], lda, &complex_zero,
                   &W[0], lda);
       memcpy(&F[0], &W[0], sizeof(complex<double>)*N*N);   
    }else if(dir==-1){
       // backward
       cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lda, lda,
                   lda, &complex_one, &F[0], lda, &eiH[0], lda, &complex_zero,
                   &W[0], lda);
       memcpy(&F[0], &W[0], sizeof(complex<double>)*N*N);   

    }else{
       printf("Error! Unknown propagation type %d in Exc::prop\n",dir);
       exit(EXIT_FAILURE);
    }
    
} 

void Exc::readHs(int nframe)
{
// read the given frame of the excitonic Hamiltonian
// from the Hamiltonian file
//
   vector<float> Htmp;
   Htmp.resize(ndim1);

   int64_t file_offset;
   
   file_offset = nframe*sizeof(float)*ndim1;
   hinfile.seekg(file_offset);

   hinfile.read(reinterpret_cast<char*>(&Htmp[0]), ndim1*sizeof(float));
   if( not hinfile.good() ) { fileReadErr(Hfile); }

   // save the one-exciton Hamiltonian:
   fill_n(H1.begin(), n1ex2, 0.0);
   int jj;
   for(int i=0; i<nchrom; ++i){
      jj = i*nchrom - i*(i+1)/2;
      for(int j=i; j<nchrom; j++){
         H1[i*nchrom + j] = static_cast<double> (Htmp[jj+j]);
         H1[j*nchrom + i] = H1[i*nchrom + j];
      }
   }

   // build two-exciton Hamiltonian for 2D IR
   if (ir2d) buildH2();

   // subtract average frequency
   for(int i=0; i<nchrom; ++i)
      H1[i*nchrom + i] -= w_avg;

   // and for the two-exciton manifold
   for(int i=0; i<n2ex; ++i)
      H2[i*n2ex + i] -= 2.0*w_avg;

}

void Exc::buildH2()
{
//
// Build two-exciton Hamiltonian for 2D IR calculations here
//
   int nind, mind;

   fill_n(H2.begin(), n2ex2, 0.0);

   // diagonal part
   for(int i=0; i<n1ex; ++i){
      nind = get2nx(i,i);
      // <i,i|H|i,i>
      H2[nind*n2ex + nind] = 2.0*H1[i*n1ex + i] - anharm;
      for(int j=i+1; j<n1ex; j++){
         nind = get2nx(i,j);
         // <i,j|H|i,j>
         H2[nind*n2ex + nind] = H1[i*n1ex + i] + H1[j*n1ex + j];
      }
   } 

   // off diagonal couplings between singly and doubly excited state
   // under harmonic approximation it is sqrt(2) of the same couplings
   // within the singly excited subspace
   for(int i=0; i<n1ex; i++){
      nind = get2nx(i,i);
      for(int j=0; j<n1ex; ++j){
         mind = get2nx(i,j);
         // <i,i|H|i,j>
         H2[nind*n2ex + mind] = sqrt(2.0)*H1[i*n1ex + j];
         H2[mind*n2ex + nind] = H2[nind*n2ex + mind];
      }
   }

   // off diagonal couplings between singly excited states <i,j|H|i,k>
   for(int i=0; i<n1ex; ++i){
      for(int j=i+1; j<n1ex; ++j){
         nind = get2nx(i,j);
         for(int k=0; k<n1ex; ++k){
            if(i==k) continue;
            mind = get2nx(i,k);
            if(nind == mind) continue;
            // <i,j|H|i,k>
            H2[nind*n2ex + mind] = H1[j*n1ex + k];
            H2[mind*n2ex + nind] = H2[nind*n2ex + mind];
         }
         for(int k=0; k<n1ex; ++k){
            if(j==k) continue;
            mind = get2nx(j,k);
            if(nind == mind) continue;
            // <i,j|H|k,j>
            H2[nind*n2ex + mind] = H1[i*n1ex + k];
            H2[mind*n2ex + nind] = H2[nind*n2ex + mind];
         }
      }
   }
   
}

int Exc::get2nx(int i, int j)
{
// return index for two-exciton states
  if(i>j){
     return nchrom*i - i*(i+1)/2 + j;
  }else{
     return nchrom*j - j*(j+1)/2 + i;
  }
}

void Exc::readDs(int nframe)
{
// read the given frame from dipole file
   int size = 3*n1ex;

   vector<float> Mtmp;
   Mtmp.resize(size);

   int64_t file_offset;
  
   file_offset = nframe*size*sizeof(float);

   //cout << file_offset << endl;
   dinfile.seekg(file_offset);

   dinfile.read(reinterpret_cast<char*>(&Mtmp[0]), size*sizeof(float));
   if( not dinfile.good() ) { fileReadErr(Dfile); }

   fill_n(mu01_x.begin(), n1ex, 0.0);
   fill_n(mu01_y.begin(), n1ex, 0.0);
   fill_n(mu01_z.begin(), n1ex, 0.0);

   for(int ii=0; ii<nchrom; ++ii){
      mu01_x[ii] = static_cast<double> (Mtmp[3*ii]);
      mu01_y[ii] = static_cast<double> (Mtmp[3*ii+1]);
      mu01_z[ii] = static_cast<double> (Mtmp[3*ii+2]);
   }

   // 2D IR
   buildM21();

}

void Exc::buildM21()
{
//
// Build two-exciton to one-exciton transition dipole moment
//
   int mx;

   fill_n(mu12_x.begin(), n1ex*n2ex, 0.0);
   fill_n(mu12_y.begin(), n1ex*n2ex, 0.0);
   fill_n(mu12_z.begin(), n1ex*n2ex, 0.0);

   for(int i=0; i<n1ex; ++i){
      mx = get2nx(i,i);
      // <i|mu|i,i>
      mu12_x[i*n2ex + mx] = sqrt(2.0)*mu01_x[i];
      mu12_y[i*n2ex + mx] = sqrt(2.0)*mu01_y[i];
      mu12_z[i*n2ex + mx] = sqrt(2.0)*mu01_z[i];
      for(int j=0; j<n1ex; ++j){
         if(j==i) continue;
         mx = get2nx(i,j);
         // <i|mu|i,j>
         mu12_x[i*n2ex + mx] = mu01_x[j];
         mu12_y[i*n2ex + mx] = mu01_y[j];
         mu12_z[i*n2ex + mx] = mu01_z[j];
      }
   } 

}

void Exc::assignD(int itime){
//
// Set transition dipole moments for a given time frame
//

   if(itime==0){
      memcpy(&mu01_t0_x[0], &mu01_x[0], sizeof(double)*n1ex);
      memcpy(&mu01_t0_y[0], &mu01_y[0], sizeof(double)*n1ex);
      memcpy(&mu01_t0_z[0], &mu01_z[0], sizeof(double)*n1ex);
   }else if(itime==1){
      memcpy(&mu01_t1_x[0], &mu01_x[0], sizeof(double)*n1ex);
      memcpy(&mu01_t1_y[0], &mu01_y[0], sizeof(double)*n1ex);
      memcpy(&mu01_t1_z[0], &mu01_z[0], sizeof(double)*n1ex);
   }else if(itime==2){
      memcpy(&mu01_t2_x[0], &mu01_x[0], sizeof(double)*n1ex);
      memcpy(&mu01_t2_y[0], &mu01_y[0], sizeof(double)*n1ex);
      memcpy(&mu01_t2_z[0], &mu01_z[0], sizeof(double)*n1ex);
      memcpy(&mu12_t2_x[0], &mu12_x[0], sizeof(double)*n1ex*n2ex);
      memcpy(&mu12_t2_y[0], &mu12_y[0], sizeof(double)*n1ex*n2ex);
      memcpy(&mu12_t2_z[0], &mu12_z[0], sizeof(double)*n1ex*n2ex);
   }else if(itime==3){
      memcpy(&mu01_t3_x[0], &mu01_x[0], sizeof(double)*n1ex);
      memcpy(&mu01_t3_y[0], &mu01_y[0], sizeof(double)*n1ex);
      memcpy(&mu01_t3_z[0], &mu01_z[0], sizeof(double)*n1ex);
      memcpy(&mu12_t3_x[0], &mu12_x[0], sizeof(double)*n1ex*n2ex);
      memcpy(&mu12_t3_y[0], &mu12_y[0], sizeof(double)*n1ex*n2ex);
      memcpy(&mu12_t3_z[0], &mu12_z[0], sizeof(double)*n1ex*n2ex);
   }else{
      printf("Error! Unknown time in Exc::assignD\n");
      exit(EXIT_FAILURE);
   }
}

void Exc::save_F1_t0t1_mu0(int it)
{
//
// Save propagated Hamiltonian and multiply by transition dipoles at t=0
//
   int i;
   MKL_INT lda;
   lda = (MKL_INT) n1ex;

   vector<complex<double>> work(n1ex,complex_zero);
   vector<complex<double>> mu(n1ex,complex_zero);

   // x
   for(i=0; i<n1ex; ++i) mu[i] = complex_one*mu01_t0_x[i];
   
   cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F1_t0t1[0], lda, &mu[0], 1, &complex_zero, &work[0], 1);
   memcpy(&F1_t0t1_mu01_0_x[it*n1ex], &work[0], sizeof(complex<double>)*n1ex);

   // y
   for(i=0; i<n1ex; ++i) mu[i] = complex_one*mu01_t0_y[i];

   cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F1_t0t1[0], lda, &mu[0], 1, &complex_zero, &work[0], 1);
   memcpy(&F1_t0t1_mu01_0_y[it*n1ex], &work[0], sizeof(complex<double>)*n1ex);

   // z
   for(i=0; i<n1ex; ++i) mu[i] = complex_one*mu01_t0_z[i];

   cblas_zgemv(CblasRowMajor, CblasNoTrans, lda, lda, &complex_one,
              &F1_t0t1[0], lda, &mu[0], 1, &complex_zero, &work[0], 1);
   memcpy(&F1_t0t1_mu01_0_z[it*n1ex], &work[0], sizeof(complex<double>)*n1ex);

}

void Exc::assignDpol(vector<complex<double>> &mu0_eg, vector<complex<double>> &mu1_eg,
                     vector<complex<double>> &mu2_eg, vector<complex<double>> &mu3_eg,
                     vector<complex<double>> &mu2_ce, vector<complex<double>> &mu3_ce,
                     int it, int k)
{
   int i;

   if(k==0){ // XXXX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }
   else if(k==1){ // YYYY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }
   else if(k==2){ //ZZZZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==3){ //XXYY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==4){ //XXZZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==5){ //YYXX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==6){ //YYZZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==7){ //ZZXX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==8){ //ZZYY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==9){ //XYXY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==10){ //XZXZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==11){ //YXYX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==12){ //YZYZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==13){ //ZXZX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==14){ //ZYZY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==15){ //XYYX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==16){ //XZZX
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_x[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_x[i];
   }else if(k==17){ //YXXY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==18){ //YZZY
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_y[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_z[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_z[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_y[i];
   }else if(k==19){ //ZXXZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_x[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_x[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_x[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else if(k==20){ //ZYYZ
      for(i=0; i<n1ex; ++i)      mu0_eg[i] = F1_t0t1_mu01_0_z[it*n1ex+i];
      for(i=0; i<n1ex; ++i)      mu1_eg[i] = mu01_t1_y[i];
      for(i=0; i<n1ex; ++i)      mu2_eg[i] = mu01_t2_y[i];
      for(i=0; i<n1ex; ++i)      mu3_eg[i] = mu01_t3_z[i];
      for(i=0; i<n1ex*n2ex; ++i) mu2_ce[i] = mu12_t2_y[i];
      for(i=0; i<n1ex*n2ex; ++i) mu3_ce[i] = mu12_t3_z[i];
   }else{
      printf("Error! Unknown value of k=%d in assignDpol.\n",k);
      exit(EXIT_FAILURE);
   }

}

complex<double> Exc::calcR2D_R1(int it)
{
   int i, k;
   complex<double> r2d, work0a, work0b;
   double weight;

   MKL_INT lda1, lda2;
   lda1 = (MKL_INT) n1ex;
   lda2 = (MKL_INT) n2ex;

   vector<complex<double>> work1a, work1b, work1c, work2a, work2b;

   work1a.resize(n1ex, complex_zero);
   work1b.resize(n1ex, complex_zero);
   work1c.resize(n1ex, complex_zero);
   work2a.resize(n2ex, complex_zero);
   work2b.resize(n2ex, complex_zero);

   r2d = complex_zero;

   // isotropically averaged spectra
   for(k=0; k<21; k++){
      assignDpol(mu0_eg, mu1_eg, mu2_eg, mu3_eg, mu2_ce, mu3_ce, it, k);

      // assign weights for parallel isotropic spectrum
      if(k <= 2) weight = 1.0/5.0;
      else weight = 1.0/15.0;

      // ===========================
      // ground state bleach:
      // ===========================
      // work1a=F1_t2t3*mu2_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &mu2_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work0a=mu3_eg*work1a
      cblas_zdotu_sub(lda1, &mu3_eg[0], 1, &work1a[0], 1, &work0a);

      // work0b=conj(mu0_eg)*mu1_eg
      for(i=0; i<n1ex; ++i) work1a[i] = conj(mu0_eg[i]);
      cblas_zdotu_sub(lda1, &work1a[0], 1, &mu1_eg[0], 1, &work0b);

      r2d += weight*work0b*work0a;

      // ===========================
      // stimulated emission:
      // ===========================
      // work1a = F1_t1t2*mu1_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu1_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work1b = F1_t2t3*work1a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &work1a[0], 1, &complex_zero, &work1b[0], 1);
      // work0a = mu3_eg*work1b
      cblas_zdotu_sub(lda1, &mu3_eg[0], 1, &work1b[0], 1, &work0a);

      // work1a = F1_t1t2*mu0_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu0_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work0b = conj(work1a)*mu2_eg
      for(i=0; i< n1ex; i++) work1b[i] = conj(work1a[i]);
      cblas_zdotu_sub(lda1, &work1b[0], 1, &mu2_eg[0], 1, &work0b);

      r2d += weight*work0b*work0a;

      // ===========================
      // excited state absorption:
      // ===========================
      // work1a = F1_t1t2*mu1_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu1_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work2a = mu2_ce'*work1a
      cblas_zgemv(CblasRowMajor, CblasTrans, lda1, lda2, &complex_one,
                  &mu2_ce[0], lda2, &work1a[0], 1, &complex_zero, &work2a[0], 1);
      // work2b = F2_t2t3*work2a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda2, lda2, &complex_one,
                  &F2_t2t3[0], lda2, &work2a[0], 1, &complex_zero, &work2b[0], 1);
      // work1c = mu3_ce*work2b
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda2, &complex_one,
                  &mu3_ce[0], lda2, &work2b[0], 1, &complex_zero, &work1c[0], 1);

      // work1a = F1_t1t2*mu0_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu0_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work1b = F1_t2t3*work1a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &work1a[0], 1, &complex_zero, &work1b[0], 1);

      // work0a = conj(work1b)*work1c
      for(i=0; i<n1ex; i++) work1a[i] = conj(work1b[i]);
      cblas_zdotu_sub(lda1, &work1a[0], 1, &work1c[0], 1, &work0a);

      r2d -= weight*work0a; 
   }

   return r2d;
}

complex<double> Exc::calcR2D_R2(int it)
{
   int i, k;
   complex<double> r2d, work0a, work0b;
   double weight;

   MKL_INT lda1, lda2;
   lda1 = (MKL_INT) n1ex;
   lda2 = (MKL_INT) n2ex;

   vector<complex<double>> work1a, work1b, work1c, work2a, work2b;

   work1a.resize(n1ex, complex_zero);
   work1b.resize(n1ex, complex_zero);
   work1c.resize(n1ex, complex_zero);
   work2a.resize(n2ex, complex_zero);
   work2b.resize(n2ex, complex_zero);


   r2d = complex_zero;

   // isotropically averaged spectra
   for(k=0; k<21; k++){

      assignDpol(mu0_eg, mu1_eg, mu2_eg, mu3_eg, mu2_ce, mu3_ce, it, k);

      // assign weights
      if(k <= 2) weight = 1.0/5.0;
         else weight = 1.0/15.0;

      // ============================
      // ground state bleach 
      // ============================
      // work1a = F1_t2t3*mu2_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &mu2_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work0a = mu3_eg*work1a
      cblas_zdotu_sub(lda1, &mu3_eg[0], 1, &work1a[0], 1, &work0a);

      // work0b = mu0_eg*mu1_eg
      cblas_zdotu_sub(lda1, &mu0_eg[0], 1, &mu1_eg[0], 1, &work0b);

      r2d += weight*work0b*work0a;

      // ============================
      // stimulated emission
      // ============================
      // work1a = F1_t1t2*mu0_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu0_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work1b = F1_t2t3*work1a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &work1a[0], 1, &complex_zero, &work1b[0], 1);
      // work0b = work1b*mu3_eg
      cblas_zdotu_sub(lda1, &work1b[0], 1, &mu3_eg[0], 1, &work0b);

      // work1a = F1_t1t3*mu1_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu1_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work0a = conj(work1a)*mu2_eg
      for(i=0; i<n1ex; ++i) work1b[i] = conj(work1a[i]);

      cblas_zdotu_sub(lda1, &work1b[0], 1, &mu2_eg[0], 1, &work0a);

      r2d += weight*work0b*work0a;

      // =============================
      // excited state absorption
      // =============================
      // work1a = F1_t1t2*mu0_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu0_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work2a = mu2_ce*work1a
      cblas_zgemv(CblasRowMajor, CblasTrans, lda1, lda2, &complex_one,
                  &mu2_ce[0], lda2, &work1a[0], 1, &complex_zero, &work2a[0], 1);
      // work2b = F2_t2t3*work2a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda2, lda2, &complex_one,
                  &F2_t2t3[0], lda2, &work2a[0], 1, &complex_zero, &work2b[0], 1);
      // work1c = mu3_ce*work2b
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda2, &complex_one,
                  &mu3_ce[0], lda2, &work2b[0], 1, &complex_zero, &work1c[0], 1);

      // work1a = F1_t1t2*mu1_eg
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t1t2[0], lda1, &mu1_eg[0], 1, &complex_zero, &work1a[0], 1);
      // work1b = F1_t2t3*work1a
      cblas_zgemv(CblasRowMajor, CblasNoTrans, lda1, lda1, &complex_one,
                  &F1_t2t3[0], lda1, &work1a[0], 1, &complex_zero, &work1b[0], 1);

      // work0a = conj(work1b)*work1c
      for(i=0; i<n1ex; ++i) work1a[i] = conj(work1b[i]);
      cblas_zdotu_sub(lda1, &work1a[0], 1, &work1c[0], 1, &work0a);

      r2d -= weight*work0a;
   }
   
   return r2d;
}

void Exc::writeR2D()
{
   string fn;
   ofstream ofile;
 
   int t1, t3;

   complex<double> Rtmp;

   // rephasing
   fn="RparI.dat";
   ofile.open(fn);
   // check if file is open

   for(t1=0; t1<nt1t3; t1++){
      for(t3=0; t3<nt1t3; t3++){
         ofile << t1*dt << " " << t3*dt << " " << R2D_R1[t1*nt1t3+t3].real() << " " << R2D_R1[t1*nt1t3+t3].imag() << endl;
      }
   }

   ofile.close();

   fn="RparII.dat";
   ofile.open(fn);

   for(t1=0; t1<nt1t3; t1++){
      for(t3=0; t3<nt1t3; t3++){
         ofile << t1*dt << " " << t3*dt << " " << R2D_R2[t1*nt1t3+t3].real() << " " << R2D_R2[t1*nt1t3+t3].imag() << endl;
      }
   }
   ofile.close();

}

void Exc::write2DRabs()
{
//
// write 2D IR spectrum
//
   string fn;
   int it1, it3, i, j;

   fftw_plan plan;

   vector<complex<double>> fftIn (NFFT*NFFT, complex_zero);
   vector<complex<double>> fftOut (NFFT*NFFT, complex_zero);
   vector<complex<double>> res (NFFT*NFFT, complex_zero);

   plan = fftw_plan_dft_2d(NFFT, NFFT, reinterpret_cast<fftw_complex*>(&fftIn[0]),
                           reinterpret_cast<fftw_complex*>(&fftOut[0]),
                           FFTW_BACKWARD, FFTW_ESTIMATE);

   printf(" done setting up plan\n");

   // Fourier transform rephasing response function
   for(i=0; i<NFFT*NFFT; i++) fftIn[i] = complex_zero;

   for(it1 = 0; it1 < nt1t3; it1++){
      for(it3 = 0; it3 < nt1t3; it3++){
         fftIn[it1*NFFT+it3] = R2D_R1[it1*nt1t3+it3];
         // divide t=0 point by 2 (see Hamm and Zanni sec. 9.5.3)
         if(it1 == 0 and it3==0) fftIn[it1*NFFT+it3] /= 2.0;
      }     
   }

   printf(" run 1st FFT \n");
   fftw_execute(plan);
   printf(" done\n");
 
   fn = "RparIw.dat";
   write2Dout(fftOut, fn, "rephasing", NFFT);
   
   // save rephasing contribution to purely absorptive spectrum
   for(i=0; i<NFFT; i++){
      for(j=0; j<NFFT; ++j){
         if(i==j) res[i*NFFT+j] = fftOut[i*NFFT+j];
         else     res[i*NFFT+j] = fftOut[(NFFT-i)*NFFT+j];
      }
   }

   // fourier transform non-rephasing response functions, see Hamm and Zanni eq 4.31
   for(i=0; i<NFFT*NFFT; ++i) fftIn[i] = complex_zero;

   for(it1 = 0; it1 < nt1t3; ++it1){
      for(it3 = 0; it3 < nt1t3; ++it3){
         fftIn[it1*NFFT+it3] = R2D_R2[it1*nt1t3+it3];
         // divide t=0 point by 2 (see Hamm and Zanni sec. 9.5.3)
         if(it1 == 0 and it3==0) fftIn[it1*NFFT+it3] /= 2.0;
      }
   }
  
   printf(" run 2nd FFT\n"); 
   fftw_execute(plan);
   printf(" done\n");

   fftw_destroy_plan(plan);

   fn = "RparIIw.dat";
   write2Dout(fftOut, fn, "non-rephasing", NFFT);

   // save rephasing contribution to purely absorptive spectrum
   for(i=0; i<NFFT; i++)
      for(j=0; j<NFFT; j++)
         res[i*NFFT+j] += fftOut[i*NFFT+j];

   fn = "RparAbs.dat";
   write2Dout(res, fn, "rabs", NFFT);

}

void Exc::write2Dout(vector<complex<double>> data, string fn, string which, int n)
{
   int i, j;
   double shift_w1, shift_w3, window0_w1, window0_w3, window1_w1, window1_w3;
   double scale, w1, w3;
   ofstream ofile;

   // scaling
   scale = dt*n/(2.0*M_PI*constants::HBAR*sqrt(n));
   scale *= -1.;
   
   ofile.open(fn);

   // assign shifts and spectral window limits
   shift_w1 = w_avg;
   shift_w3 = w_avg;
   window0_w1 = 1400.0;
   window1_w1 = 1700.0;
   window0_w3 = 1400.0;
   window1_w3 = 1700.0;

   if(which.compare("rephasing")==0){
      printf(" Rephasing parallel ZZZZ polarized spectrum\n");
      shift_w1 = -w_avg;
      window0_w1 = -1400.0; //window1;
      window1_w1 = -1700.0; //window0;
   }
   else if(which.compare("non-rephasing")==0){
      printf(" Non-rephasing parallel ZZZZ polarized spectrum\n");
   }
   else if(which.compare("rabs")==0){
      printf(" Pure parallel ZZZZ polarized spectrum\n");
   }
   else{
      printf(" Error in write2Dout\n");
      exit(EXIT_FAILURE);
   }

   for(i = n/2; i<n; i++){
      w1 = 2.0*M_PI*constants::HBAR*(i-n)/(n*dt) + shift_w3;
      if(w1 < window0_w1 or w1 > window1_w1) continue;
      for(j = n/2; j<n; j++){
         w3 = 2.0*M_PI*constants::HBAR*(j-n)/(dt*n) + shift_w3;
         if(w3 < window0_w3 or w3 > window1_w3) continue;
         ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() << " " << scale*data[i*n+j].imag() << endl;
      }
      for(j=0; j<n/2; ++j){
         w3 = 2.0*M_PI*constants::HBAR*j/(dt*n) + shift_w3;
         if(w3 < window0_w3 or w3 > window1_w3) continue;
         ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() << " " << scale*data[i*n+j].imag() << endl;
      }
   }

   for(i =0; i<n/2; i++){
      w1 = 2.0*M_PI*constants::HBAR*i/(n*dt) + shift_w1;
      if(w1 < window0_w1 or w1 > window1_w1) continue;
      for(j = n/2; j<n; j++){
         w3 = 2.0*M_PI*constants::HBAR*(j-n)/(dt*n) + shift_w3;
         if(w3 < window0_w3 or w3 > window1_w3) continue;
         ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() << " " << scale*data[i*n+j].imag() << endl;
      }
      for(j=0; j<n/2; ++j){
         w3 = 2.0*M_PI*constants::HBAR*j/(dt*n) + shift_w3;
         if(w3 < window0_w3 or w3 > window1_w3) continue;
         ofile << w1 << " " << w3 << " " << scale*data[i*n+j].real() << " " << scale*data[i*n+j].imag() << endl;
      }
   }
   ofile.close();


}  
