#include "wmaps.h"

using namespace std;

WmapS::WmapS(string map_name_inp, string chrom_type, bool intramolOHc) :
             map_name_inp(map_name_inp), ctype(chrom_type), intrac(intramolOHc)
{
   printf("\n** Reading Stretch Map **\n");
   map_name = str_toupper(map_name_inp);

   if(map_name=="LI_2010_TIP4P"){
      printf("   Using hydroxyl stretch map from : F. Li and J. L. Skinner, J. Chem. Phys. 132, 244504 (2010)\n");
      printf("   Note! This map can only describe 0->1 transitions.\n");
      emap = {3732.90, -3519.8, -153520.0, 
              0.19318, -1.7248e-5, 0.0, 
              1.6102, 5.8697e-4, 0.0, 
              0.1622, 10.381, 137.6, 
              2748.2, -2572.2, -102980.0,
              0.16598, -2.0752e-5, 0.0,
              1.9813, 9.1419e-4, 0.0,
              0.1622, 10.381, 137.6,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
             -1361.0, 27165.0, -1.887};
   }else if(map_name=="GRUENBAUM_2013_TIP4P"){
      printf("   Using hydroxyl stretch map from : S. M. Gruenbaum et al., J. Chem. Theory Comput. 9, 3109 (2013)\n");
      emap = {3760.2,  -3541.7,     -152677.0,
              0.19285, -1.7261e-5,   0.0,
              1.6466,   5.7692e-4,   0.0,
              0.1646,   11.39,       63.41,
              2767.8,  -2630.3,     -102601.0,
              0.16593, -2.0632e-5,   0.0,
              2.0475,   8.9108e-4,   0.0,
              0.1646,   11.39,       63.41,
              3606.0,  -3498.6,     -198715.0,
              0.26836, -2.3788e-5,   0.0,
              2.0160,   8.7684e-4,   0.0,
              0.1646,   11.39,       63.41,
              2673.0,  -1763.5,     -138534.0,
              0.23167, -2.8596e-5,   0.0,
              2.6233,   13.1443e-4,  0.0,
              0.1646,   11.39,       63.41,
             -1361.0,   27165.0,    -1.887}; 
   }else if(map_name=="AUER_2008_SPCE"){
      printf("   Using hydroxyl stretch map from : B. M. Auer et al., J. Chem. Phys. 128, 224511 (2008)\n");
      emap = {3762.0,  -5060.0,    -86225.0,
              0.1934,  -1.75e-5,    0.0,
              1.611,    5.893e-4,   0.0,
              0.1333,   14.17,      0.0,
              2762.6,  -3640.8,    -56641.0,
              0.16627, -2.0884e-5,  0.0,
              1.9844,   9.1907e-4,  0.0,
              0.1333,   14.17,      0.0,
              3614.1,  -5493.7,    -115670.0,
              0.1428,  -1.29e-5,    0.0,
              0.0,      0.0,        0.0,
              0.1333,   14.17,      0.0,
              2695.8,  -3785.1,    -73074.0,
              0.1229,  -1.525e-5,   0.0,
              0.0,      0.0,        0.0,
              0.1333,   14.17,      0.0,
             -1789.0,   23852.0,   -1.966}; 
   }else{
      printf(" Error! Cannot recognize the map! %s \n",map_name.c_str());
      exit(EXIT_FAILURE);
   };

   if(!intrac)
      printf("   Note! OH-stretch intramolecular couplings are set to 0.\n");

}

float WmapS::getw01E(float E){ 
   float rv=0.0;
   if(ctype=="wsOH" || ctype=="wswbH2O"){
      rv = getw01E_OH(E);
   }else if(ctype=="wsOD" || ctype=="wswbD2O"){
      rv = getw01E_OD(E);
   }
   return rv;
}

float WmapS::getm01E(float E){
   float rv = 0.0;
   if(ctype=="wsOH" || ctype=="wswbH2O"){
      rv = getm01E_OH(E);
   }else if(ctype=="wsOD" || ctype=="wswbD2O"){
      rv = getm01E_OD(E);
   }
   return rv;
}

float WmapS::getx01E(float E){
   float rv = 0.0;
   if(ctype=="wsOH" || ctype=="wswbH2O"){
      rv = getx01E_OH(E);
   }else if(ctype=="wsOD" || ctype=="wswbD2O"){
      rv = getx01E_OD(E);
   }
   return rv;
}

float WmapS::getp01E(float E){
   float rv = 0.0;
   if(ctype=="wsOH" || ctype=="wswbH2O"){
      rv = getp01E_OH(E);
   }else if(ctype=="wsOD" || ctype=="wswbD2O"){
      rv = getp01E_OD(E);
   }
   return rv;
}

float WmapS::getw01E_OH(float E)
{ return emap.w10_OH_0 + emap.w10_OH_1*E + emap.w10_OH_2*E*E; }

float WmapS::getw01E_OD(const float E)
{ return emap.w10_OD_0 + emap.w10_OD_1*E + emap.w10_OD_2*E*E; }

float WmapS::getm01E_OH(const float E)
{ return emap.m10_OH_0 + emap.m10_OH_1*E + emap.m10_OH_2*E*E; }

float WmapS::getm01E_OD(const float E)
{ return emap.m10_OD_0 + emap.m10_OD_1*E + emap.m10_OD_2*E*E; }

float WmapS::getp01E_OH(const float w)
{ return emap.p10_OH_0 + emap.p10_OH_1*w; }

float WmapS::getp01E_OD(const float w)
{ return emap.p10_OD_0 + emap.p10_OD_1*w; }

float WmapS::getx01E_OH(const float w)
{ return emap.x10_OH_0 + emap.x10_OH_1*w; }

float WmapS::getx01E_OD(const float w)
{ return emap.x10_OD_0 + emap.x10_OD_1*w; }

float WmapS::getcnn(const float Ei, const float Ej)
{
//
//  Intramolecular coupling of two hydroxyl stretches
// 
   float xi = 0.0; 
   float xj = 0.0;
   float pi = 0.0;
   float pj = 0.0;
   float wi = 0.0;
   float wj = 0.0;
   float wc = 0.0;

   if(intrac){
      if(ctype=="wsOH"){
         wi = getw01E_OH(Ei);
         wj = getw01E_OH(Ej);
         xi = getx01E_OH(wi);
         xj = getx01E_OH(wj);
         pi = getp01E_OH(wi);
         pj = getp01E_OH(wj);
      }else if(ctype=="wsOD"){
         wi = getw01E_OD(Ei);
         wj = getw01E_OD(Ej);
         xi = getx01E_OD(wi);
         xj = getx01E_OD(wj);
         pi = getp01E_OD(wi);
         pj = getp01E_OD(wj);
      }
      wc = (emap.wij0 + emap.wije*(Ei + Ej))*xi*xj + emap.wijpp*pi*pj;
   }
   return wc;
}

float WmapS::getcnn(const float ei, const float xi, const float pi,
                    const float ej, const float xj, const float pj)
{
   double wc;
   wc = (emap.wij0 + emap.wije*(ei + ej))*xi*xj + emap.wijpp*pi*pj;
   return wc;
} 
