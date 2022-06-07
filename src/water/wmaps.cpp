#include "wmaps.h"

using namespace std;

Wmap::Wmap(string map_name)
{
    printf("\n** Reading Spectroscopic Map **\n");

    if(map_name=="skinner_tip4p_2010_OH"){
       printf("   Using OH stretch map from : F. Li and J. L. Skinner, J. Chem. Phys. 132, 244504 (2010)\n");
       emap = {3732.90, -3519.8, -1.5352e5, 
               0.19318, -1.7248e-5, 0.0, 
               1.6102, 5.8697e-4, 0.0, 
               0.1622, 10.381, 137.6, 
              -1361.0, 27165.0, -1.887};
    }else if(map_name=="skinner_tip4p_2010_OD"){
       printf("   Using OD stretch map from : F. Li and J. L. Skinner, J. Chem. Phys. 132, 244504 (2010)\n");
       emap = {2748.2, -2572.2, -102980.0,
               0.16598, -2.0752e-5, 0.0,
               1.9813, 9.1419e-4, 0.0,
               0.1622, 10.381, 137.6,
              -1361.0, 27165.0, -1.887};
    }else if(map_name=="gruenbaum_tip4p_2013_OH"){
       printf("   Using OH stretch map from : S. M. Gruenbaum et al., J. Chem. Theory Comput. 9, 3109 (2013)\n");
       emap = {3760.2,  -3541.7, -152677,
               0.19285, -1.7261e-5, 0.0,
               1.6466,   5.7692e-4, 0.0,
               0.1646,  11.39, 63.41,
               -1361.0, 27165.0, -1.887}; 
    }else if(map_name=="gruenbaum_tip4p_2013_OD"){
       printf("   Using OD stretch map from : S. M. Gruenbaum et al., J. Chem. Theory Comput. 9, 3109 (2013)\n");
       emap = {2767.8, -2630.3, -102601,
               0.16593, -2.0632e-5, 0.0,
               2.0475, 8.9108e-4, 0.0,
               0.1646,  11.39, 63.41,
               -1361.0, 27165.0, -1.887}; 
   };
    

}

float Wmap::getw01E(const float E)
{ return emap.w0 + emap.w1*E + emap.w2*E*E; }

float Wmap::getm01E(const float E)
{ return emap.m0 + emap.m1*E + emap.m2*E*E; }

float Wmap::getp01E(const float E)
{ return emap.p0 + emap.p1*E + emap.p2*E*E; }

float Wmap::getx01E(const float E)
{ return emap.x0 + emap.x1*E + emap.x2*E*E; }

float Wmap::getcnn(const float Ei, const float Ej)
{
   float xi, xj, pi, pj, w, wi, wj;
   wi = getw01E(Ei);
   wj = getw01E(Ej);
   xi = getx01E(wi);
   xj = getx01E(wj);
   pi = getp01E(wi);
   pj = getp01E(wj);
   w = (emap.wij0 + emap.wije*(Ei + Ej))*xi*xj + emap.wijpp*pi*pj;
   return w;
}
