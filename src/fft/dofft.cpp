#include "dofft.h"

using namespace std;

void FFT1D(vector<complex<double>> fft_in, vector<complex<double>> &fft_out, 
           double dt, int nfft)
{

   complex<double> complex_zero = {0.0, 0.0};
   complex<double> img = {0.0, 1.0};

   double scale;
   int ntime;

   fftw_plan plan;
   
   fill_n(fft_out.begin(), nfft, complex_zero);
   vector<complex<double>> fftIn (nfft, complex_zero);
   vector<complex<double>> fftOut(nfft, complex_zero);

   scale = dt*nfft/(2.0*M_PI*constants::HBAR*sqrt(nfft));
  
   plan = fftw_plan_dft_1d(nfft, reinterpret_cast<fftw_complex*>(&fftIn[0]),
                           reinterpret_cast<fftw_complex*>(&fftOut[0]),
                           FFTW_FORWARD, FFTW_MEASURE);

   ntime = fft_in.size();

   if(2*ntime > nfft){
      printf(" Increase NFFT to at least %d ",2*ntime+1);
      exit(EXIT_FAILURE);
   }

   for(int it=0; it<ntime; ++it){
      fftIn[it] = fft_in[it];
      if(it==0) continue;
      fftIn[nfft-it] = conj(fft_in[it]);  
   }

   fftw_execute(plan);

   for(int it=0; it<nfft; ++it)
      fft_out[it].real(fftOut[it].real());

   fill_n(fftIn.begin(), nfft, complex_zero);
   for(int it=0; it<ntime; ++it){
      fftIn[it] = img*fft_in[it];
      if(it==0) continue;
      fftIn[nfft-it] = conj(img*fft_in[it]);
   } 

   fftw_execute(plan);

   for(int it=0; it<nfft; ++it)
      fft_out[it].imag(fftOut[it].real());

   for(int it=0; it<nfft; ++it)
      fft_out[it] *= scale;

   fftw_destroy_plan(plan);
}

