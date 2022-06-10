#include "wmaps.h"

using namespace std;

Wmap::Wmap(string map_name, string chrom_type) :
           ctype(chrom_type)
{
   printf("\n** Reading Spectroscopic Map **\n");

   if(map_name=="skinner_tip4p_2010"){
      printf("   Using hydroxyl stretch map from : F. Li and J. L. Skinner, J. Chem. Phys. 132, 244504 (2010)\n");
      emap = {3732.90, -3519.8, -1.5352e5, 
              0.19318, -1.7248e-5, 0.0, 
              1.6102, 5.8697e-4, 0.0, 
              0.1622, 10.381, 137.6, 
              2748.2, -2572.2, -102980.0,
              0.16598, -2.0752e-5, 0.0,
              1.9813, 9.1419e-4, 0.0,
              0.1622, 10.381, 137.6,
             -1361.0, 27165.0, -1.887};
   }else if(map_name=="gruenbaum_tip4p_2013"){
      printf("   Using hydroxyl stretch map from : S. M. Gruenbaum et al., J. Chem. Theory Comput. 9, 3109 (2013)\n");
      emap = {3760.2,  -3541.7, -152677,
              0.19285, -1.7261e-5, 0.0,
              1.6466,   5.7692e-4, 0.0,
              0.1646,  11.39, 63.41,
              2767.8, -2630.3, -102601,
              0.16593, -2.0632e-5, 0.0,
              2.0475, 8.9108e-4, 0.0,
              0.1646,  11.39, 63.41,
             -1361.0, 27165.0, -1.887}; 
   }else{
      printf(" Error! Cannot recognize the map! %s \n",map_name.c_str());
      exit(EXIT_FAILURE);
   };

   // perform some checks
   if(ctype!="wsOH")
     if(ctype!="wsOD")
       if(ctype!="wsiso"){ 
          printf(" Error! Cannot recognize spec_type %s \n",ctype.c_str());
          exit(EXIT_FAILURE);
       }

}

float Wmap::getw01E(float E){ 
   if(ctype=="wsOH")
   {
      return getw01E_OH(E);
   }
   else if(ctype=="wsOD")
   {
      return getw01E_OD(E);
   }
}

float Wmap::getm01E(float E){
   if(ctype=="wsOH")
   {
      return getm01E_OH(E);
   }
   else if(ctype=="wsOD")
   {
      return getm01E_OD(E);
   }
}

float Wmap::getx01E(float E){
   if(ctype=="wsOH")
   {
      return getx01E_OH(E);
   }
   else if(ctype=="wsOD")
   {
      return getx01E_OD(E);
   }
}

float Wmap::getp01E(float E){
   if(ctype=="wsOH")
   {
      return getp01E_OH(E);
   }
   else if(ctype=="wsOD")
   {
      return getp01E_OD(E);
   }
}

float Wmap::getw01E_OH(float E)
{ return emap.w0_OH + emap.w1_OH*E + emap.w2_OH*E*E; }

float Wmap::getw01E_OD(const float E)
{ return emap.w0_OD + emap.w1_OD*E + emap.w2_OD*E*E; }

float Wmap::getm01E_OH(const float E)
{ return emap.m0_OH + emap.m1_OH*E + emap.m2_OH*E*E; }

float Wmap::getm01E_OD(const float E)
{ return emap.m0_OD + emap.m1_OD*E + emap.m2_OD*E*E; }

float Wmap::getp01E_OH(const float E)
{ return emap.p0_OH + emap.p1_OH*E + emap.p2_OH*E*E; }

float Wmap::getp01E_OD(const float E)
{ return emap.p0_OD + emap.p1_OD*E + emap.p2_OD*E*E; }

float Wmap::getx01E_OH(const float E)
{ return emap.x0_OH + emap.x1_OH*E + emap.x2_OH*E*E; }

float Wmap::getx01E_OD(const float E)
{ return emap.x0_OD + emap.x1_OD*E + emap.x2_OD*E*E; }

float Wmap::getcnn(const float Ei, const float Ej)
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

   if(ctype=="wsOH"){
      wi = getw01E_OH(Ei);
      wj = getw01E_OH(Ej);
      xi = getx01E_OH(wi);
      xj = getx01E_OH(wj);
      pi = getp01E_OH(wi);
      pj = getp01E_OH(wj);
   }
   else if(ctype=="wsOD")
   {
      wi = getw01E_OD(Ei);
      wj = getw01E_OD(Ej);
      xi = getx01E_OD(wi);
      xj = getx01E_OD(wj);
      pi = getp01E_OD(wi);
      pj = getp01E_OD(wj);
   }
   wc = (emap.wij0 + emap.wije*(Ei + Ej))*xi*xj + emap.wijpp*pi*pj;
   return wc;
}
