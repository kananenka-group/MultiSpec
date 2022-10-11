#include "nncmap.h"

using namespace std;

NNCmap::NNCmap(string map_name_inp) :
             map_name_inp(map_name_inp)
{
   map_name = str_toupper(map_name_inp);
   printf("\n** Reading nearest-neighbor map. **\n");

   if(map_name=="JANSEN_2006"){
       printf("   Using nearest-neighbor coupling map from : Jansen et al. J. Chem. Phys. 125, 044312 (2006).\n");
       // see Supporting info for the paper above.
       coupling = { 6.331, 5.951, 4.853, 5.668, 9.597,13.966,21.427,13.966, 9.597, 5.668, 4.853, 5.951, 6.331, 
                    5.645, 4.666, 5.023, 8.648,13.250,18.309,15.142, 9.130, 4.377, 2.899, 4.712, 6.081, 5.645, 
                    4.333, 3.457, 5.621,10.432,15.718,18.485, 9.208, 2.179,-1.699, 0.142, 3.562, 4.969, 4.333, 
                    2.604, 2.194, 5.158,10.530,15.330,11.107, 2.546,-5.915,-6.721,-2.935, 1.790, 3.453, 2.604, 
                    0.248, 0.513, 3.309, 8.556, 9.372, 3.675,-7.395,-11.367,-8.648,-4.273, 0.026, 0.904, 0.248,
                   -1.219,-0.279, 2.404, 4.462, 2.328,-5.560,-12.700,-12.308,-7.319,-2.714,-0.361,-0.648,-1.219,
                   -1.948,-0.646, 0.771, 0.653,-3.904,-10.815,-14.526,-10.815,-3.904, 0.653, 0.771,-0.646,-1.948,
                   -1.219,-0.648,-0.361,-2.714,-7.319,-12.308,-12.700,-5.560, 2.328, 4.462, 2.404,-0.279,-1.219, 
                    0.248, 0.904, 0.026,-4.273,-8.648,-11.367,-7.395, 3.675, 9.372, 8.556, 3.309, 0.513, 0.248, 
                    2.604, 3.453, 1.790,-2.935,-6.721,-5.915, 2.546,11.107,15.330,10.530, 5.158, 2.194, 2.604, 
                    4.333, 4.969, 3.562, 0.142,-1.699, 2.179, 9.208,18.485,15.718,10.432, 5.621, 3.457, 4.333, 
                    5.645, 6.081, 4.712, 2.899, 4.377, 9.130,15.142,18.309,13.250, 8.648, 5.023, 4.666, 5.645, 
                    6.331, 5.951, 4.853, 5.668, 9.597,13.966,21.427,13.966, 9.597, 5.668, 4.853, 5.951, 6.331};
       nTheta = (int) sqrt(coupling.size());
       dTheta = 360.0/(nTheta-1);
   }else{
      printf("Error! No NN amide I map coupling is specified. Available maps: 'Jansen_2006'.\n");      
      exit(EXIT_FAILURE);
   }

}

float NNCmap::getCoupling(const float &x, const float &y){
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
  z1=coupling[ny    *nTheta + nx  ];
  z2=coupling[(ny+1)*nTheta + nx  ];
  z3=coupling[(ny+1)*nTheta + nx+1];
  z4=coupling[ny    *nTheta + nx+1];
  
  float u,t;
  u=(x-xl)/(xh-xl);
  t=(y-yl)/(yh-yl);
  
  float ans=(1-t)*(1-u)*z1 + t*(1-u)*z2 + t*u*z3 + (1-t)*u*z4;
  return ans;
}

