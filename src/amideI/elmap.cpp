#include "elmap.h"

using namespace std;

ELmap::ELmap(string map_name_inp, string res_map_file) :
             map_name_inp(map_name_inp), res_map_file(res_map_file)
{
   map_name = str_toupper(map_name_inp);
   printf("\n** Reading backbone map. **\n");

   if(map_name=="WANG_2011"){
       printf("   Using backbone map from : Wang et al. J. Phys. Chem. B 115, 3713-3724 (2011).\n");
       mapi = 1;
   }else if(map_name=="LIN_2009"){
       printf("   Using backbone map from : Lin et al. J. Phys. Chem. B 113, 592-602 (2009).\n");
       mapi = 2;
   }else{
      printf("Error! No backbone amide I map is specified. Available maps: 'Wang_2011'.\n");      
      exit(EXIT_FAILURE);
   }

}

float ELmap::getw01(const float Ec, const float En){
  float w01 = 0.0;
  if(mapi==1){
     w01 = 1684.0 + 7729.0*Ec - 3576.0*En;
  }else if(mapi==2){
     w01 = 1717.0 + 4213.0*Ec + 2108.0*En;
  }
  return w01;
}

void ELmap::checkResMaps()
{
   if(!res_map_file.empty()){

      ifstream file(res_map_file);
      if(!file.good()){
         printf("Error! cannot open residues_maps_file %s \n",res_map_file.c_str());
         exit(EXIT_FAILURE);
      }

   
      //printf("   Reading residues and corresponding maps.\n");

      //string line, entry, substr;
 
      //int counter=0;

      //while(getline(file, line)) {
      //   counter++;
      //   stringstream linestream(line);
      // not finished...

      file.close();
   }
   
}

