#include "nummath.h"

using namespace std;

double simpsonInt(vector<double> x, vector<complex<double>> y)
{

    int N = x.size() - 1;
    double h = x[1] - x[0]; 

    double result = 0.0;

    for(int i=1; i<N; i=i+2)
       result += (h/3.0)*(y[i-1].real() + 4.0*y[i].real() + y[i+1].real()); 
    
    if(N % 2 == 1){
        result += y[N].real()*5.0*h/12.0;
        result += y[N-1].real()*2.0*h/3.0;
        result -= y[N-2].real()*h/12.0;
    }
    return result;
}
