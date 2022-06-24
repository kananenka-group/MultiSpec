#include "wmapb.h"

using namespace std;

WmapB::WmapB(string map_name, string chrom_type) :
             ctype(chrom_type)
{
   if((ctype=="wswbH2O") || (ctype=="wswbD2O") || (ctype=="wswbiso")){
      printf("\n** Reading Bend Map **\n");
      if(map_name=="ni_2015_tip4p"){
         printf("   Using bending map from : Y. Ni and J. L. Skinner, J. Chem. Phys. 143, 014502 (2015)\n");
         emap = {1581.46, 2938.51,
                 1551.32, 3147.80,
                 0.220276, -4.23217e-5,
                 0.31143, 6.00186e-5,
                -0.32163,
                -0.220589, 0.544310, -0.871378,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0,
                 0.0, 0.0, 0.0,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0, 0.0,
                 0.0,
                 0.0, 0.0, 0.0};
      }else if(map_name=="ni_2015_kananenka_2019_tip4p"){
         printf("   Using bending map from : Y. Ni and J. L. Skinner, J. Chem. Phys. 143, 014502 (2015)\n");
         printf("   updated with D2O and HOD bending frequencies taken from J. Phys. Chem. B 117, 15319 (2013)\n");
         printf("   and used in J. Phys. Chem. B 123, 5139-5146 (2019).\n");
         emap = {1581.46, 2938.51,
                 1551.32, 3147.80,
                 0.220276, -4.23217e-5,
                 0.31143, 6.00186e-5,
                -0.32163,
                -0.220589, 0.544310, -0.871378,
                 1210.0, 0.0,      // field-independent fundamental bend map
                 1170, 0.0,        // w21 with 40 cm-1 anharmonicity 
                 0.0, 0.0,
                 0.0, 0.0,
                -0.32163,
                -0.220589, 0.544310, -0.871378,
                 1465.0, 0.0,      // HOD bend fundamental see JPCB 117, 15319 (2013), Fig. 2
                 1485.0, 0.0,      // HOD 1->2 taken such that HOD overtones matches JPCB 117, 15319 (2013)
                 0.0, 0.0,
                 0.0, 0.0,
                -0.32163,
                -0.220589, 0.544310, -0.871378};
      }else{
         printf(" Error! Cannot recognize the map! %s \n",map_name.c_str());
         exit(EXIT_FAILURE);
      }
   }

}

float WmapB::getw02E(float E){  
   if(ctype=="wswbH2O"){ 
     return getw01E_HOH(E) + getw12E_HOH(E);
   }else if(ctype=="wswbD2O"){
     return getw01E_DOD(E) + getw12E_DOD(E);
   }
}

float WmapB::getw01E_HOH(float E)
{ return emap.w0_HOH_10 + emap.w1_HOH_10*E; }

float WmapB::getw01E_DOD(float E)
{ return emap.w0_DOD_10 + emap.w1_DOD_10*E; }

float WmapB::getw12E_HOH(float E)
{ return emap.w0_HOH_21 + emap.w1_HOH_21*E; }

float WmapB::getw12E_DOD(float E)
{ return emap.w0_DOD_21 + emap.w1_DOD_21*E; }

float WmapB::getw01E_HOD(float E)
{ return emap.w0_HOD_10 + emap.w1_HOD_10*E; }

float WmapB::getw12E_HOD(float E)
{ return emap.w0_HOD_21 + emap.w1_HOD_21*E; }

