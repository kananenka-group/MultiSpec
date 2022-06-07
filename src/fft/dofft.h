#ifndef FFT_H
#define FFT_H

#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include "../const/const.h"

using namespace std;

void FFT1D(vector<complex<double>> fft_in, vector<complex<double>> &fft_out, 
           double dt, int nfft);

#endif
