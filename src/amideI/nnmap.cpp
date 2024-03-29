#include "nnmap.h"

using namespace std;

NNmap::NNmap(string map_name_inp) :
             map_name_inp(map_name_inp)
{
   map_name = str_toupper(map_name_inp);
   printf("\n** Reading nearest-neighbor map. **\n");

   if(map_name=="JANSEN_2006"){
       printf("   Using nearest-neighbor map from : Jansen et al. J. Chem. Phys. 125, 044312 (2006).\n");
       // NN frequency shifts, copied from Jansen's AmideImap program
       // see Supporting info for the paper above.
       NtermShift ={ 0.404,-1.205, 0.789,13.541,22.541,28.548,50.551,28.548,22.541,13.541, 0.789,-1.205, 0.404,
                    -1.887,-1.141, 5.762,14.643,24.000,34.820,24.308,17.676,13.092, 4.738,-2.997,-2.603,-1.887,
                    -3.387,-0.124, 9.386,14.208,22.015,30.925,14.221,11.718, 4.889, 1.667,-3.653,-4.087,-3.387,
                    -0.069, 8.421,13.625,15.852,21.411,21.014,14.828,-0.075,-0.872, 6.685, 0.784,-1.270,-0.069,
                     6.132,14.360,15.442,17.498,21.387,17.706,-2.303,-8.094,-4.944, 4.026,-7.002, 2.536, 6.132,
                    10.251, 3.635, 9.426,13.395,13.176,-0.075,-7.130,-7.167,-3.443, 4.315, 3.121, 5.735,10.251,
                    -7.081, 7.816, 6.365, 5.061, 1.329,-4.588,-5.550,-4.588, 1.329, 5.061, 6.365, 7.816,-7.081,
                    10.251, 5.735, 3.121, 4.315,-3.443,-7.167,-7.130,-0.075,13.176,13.395, 9.426, 3.635,10.251, 
                     6.132, 2.536,-7.002, 4.026,-4.944,-8.094,-2.303,17.706,21.387,17.498,15.442,14.360, 6.132,
                    -0.069,-1.270, 0.784, 6.685,-0.872,-0.075,14.828,21.014,21.411,15.852,13.625, 8.421,-0.069,
                    -3.387,-4.087,-3.653, 1.667, 4.889,11.718,14.221,30.925,22.015,14.208, 9.386,-0.124,-3.387,
                    -1.887,-2.603,-2.997, 4.738,13.092,17.676,24.308,34.820,24.000,14.643, 5.762,-1.141,-1.887,
                     0.404,-1.205, 0.789,13.541,22.541,28.548,50.551,28.548,22.541,13.541, 0.789,-1.205, 0.404};
       CtermShift ={-13.410,-9.262,-0.184, 8.535,17.984,30.352,48.440,30.352,17.984, 8.535,-0.184,-9.262,-13.410,
                    -10.729,-3.802, 3.846,14.153,26.874,41.926,34.992,16.279, 5.857, 1.667,-5.281,-11.611,-10.729,
                     -4.913, 0.205, 4.514,16.736,31.155,35.848,19.171, 5.586,-5.799,-8.260,-8.499,-8.982,-4.913,
                     -0.695,-0.631, 6.910,16.212,24.955,20.718,10.336,-4.799,-20.484,-22.949,-10.223,-4.079,-0.695,
                      3.455, 2.482, 7.735,15.208,13.495, 9.584,-4.401,-22.259,-30.507,-27.853, 2.842, 1.872, 3.455,
                     10.288,16.772,16.805,13.984, 7.758,-6.108,-18.810,-25.888,-28.228,-13.221, 5.214, 9.661,10.288,
                      1.788,15.052,12.962, 6.891,-11.364,-22.800,-22.201,-22.800,-11.364, 6.891,12.962,15.052, 1.788,
                     10.288, 9.661, 5.214,-13.221,-28.228,-25.888,-18.810,-6.108, 7.758,13.984,16.805,16.772,10.288,
                      3.455, 1.872, 2.842,-27.853,-30.507,-22.259,-4.401, 9.584,13.495,15.208, 7.735, 2.482, 3.455,
                     -0.695,-4.079,-10.223,-22.949,-20.484,-4.799,10.336,20.718,24.955,16.212, 6.910,-0.631,-0.695,
                     -4.913,-8.982,-8.499,-8.260,-5.799, 5.586,19.171,35.848,31.155,16.736, 4.514, 0.205,-4.913,
                    -10.729,-11.611,-5.281, 1.667, 5.857,16.279,34.992,41.926,26.874,14.153, 3.846,-3.802,-10.729,
                    -13.410,-9.262,-0.184, 8.535,17.984,30.352,48.440,30.352,17.984, 8.535,-0.184,-9.262,-13.410};
       nTheta = (int) sqrt(NtermShift.size());
       dTheta = 360.0/(nTheta-1);
   }else{
      printf("Error! No NN amide I map is specified. Available maps: 'Jansen_2006'.\n");      
      exit(EXIT_FAILURE);
   }

}

float NNmap::getNNshift(const float &x, const float &y, const string termshift){
//
// modified Steve Strong's code
// 
  int nx=(int) (x+180.0)/dTheta;
  int ny=(int) (y+180.0)/dTheta;

  //account for case when x or y == +180
  if (nx==nTheta-1) nx--;
  if (ny==nTheta-1) ny--;
  
  if (nx<0 || nx>=nTheta || ny<0 || ny>=nTheta) {
     printf("Error! Ramachandran angles out of bounds.\n");
     exit(EXIT_FAILURE);
  }
  
  //corner points surrounding query point
  float xl,xh,yl,yh;
  xl=nx*dTheta-180.0;
  xh=xl+dTheta;
  yl=ny*dTheta-180.0;
  yh=yl+dTheta;
  
  //shift values at corners
  float z1=0;
  float z2=0;
  float z3=0;
  float z4=0;
  if(termshift.compare("Cterm")==0)
  {
     z1=CtermShift[ny    *nTheta + nx  ];
     z2=CtermShift[(ny+1)*nTheta + nx  ];
     z3=CtermShift[(ny+1)*nTheta + nx+1];
     z4=CtermShift[ny    *nTheta + nx+1];
  }
  else if(termshift.compare("Nterm")==0)
  {
     z1=NtermShift[ny    *nTheta + nx  ];
     z2=NtermShift[(ny+1)*nTheta + nx  ];
     z3=NtermShift[(ny+1)*nTheta + nx+1];
     z4=NtermShift[ny    *nTheta + nx+1];
  }else{
     printf("Error in getNNshift.\n");
     exit(EXIT_FAILURE);
  }
  
  float u,t;
  u=(x-xl)/(xh-xl);
  t=(y-yl)/(yh-yl);
  
  float ans=(1-t)*(1-u)*z1 + t*(1-u)*z2 + t*u*z3 + (1-t)*u*z4;
  return ans;
}

