#include "amideI.h"

amideI::amideI(string gro_file, string traj_file, vector<string> itp_files,
               string top_file, string spec_type, vector<string> isolabels) :
               s(gro_file, itp_files, top_file), traj_file(traj_file), 
               jobType(spec_type), isolabels(isolabels)
{
//
// At this point all molecular information has been processed
//
   atoms = s.getSystemData();
   natoms= s.getNatoms();   

   findAmideI();
   amideIJob();


}

amideI::~amideI() { }

void amideI::findAmideI()
{
//
// Here we will create a list of all amide I groups in the system
// assuming the following order:
//
   printf("\n** Searching for peptide groups assuming 'C' 'O' 'N' 'H' order. **\n");
   namideI = 0;
   for(int i=0; i<(natoms-3); ++i){
     if(atoms[i].atomName.compare("C")==0   &&
        atoms[i+1].atomName.compare("O")==0 &&
        atoms[i+2].atomName.compare("N")==0 &&
        atoms[i+3].atomName.compare("H")==0){
           namideI++;
           amideI_Clist.push_back(i);
        }
   }

   if(namideI>0){
      printf("   Found %d peptide groups.\n",namideI);
   }else{
      printf("Error! Did not find any peptide groups.\n");
      exit(EXIT_FAILURE);
   }

}

void amideI::amideIJob(){
//
// Decide what type of job was requested by the user
//
   int nchroma;
   printf("\n** Reading job type **\n");
   if(jobType=="full"){
      printf("   All peptide groups will be used to calculate spectra.\n");
      nchrom = namideI;
      chrom_Clist.assign(amideI_Clist.begin(), amideI_Clist.end());
   }else if(jobType=="iso"){
      printf("   Only isotope labeled residues will contribute to the spectra.\n");
      printf("   The following residues are isotope labeled: ");
      for(unsigned int i=0;  i<isolabels.size();++i)
         printf("%s ",isolabels[i].c_str());
      printf("\n");
      // loop over all atoms creating a chromophoe list 
      printf("   Selecting C=O groups: \n");
      for(unsigned int j=0; j<isolabels.size();++j){
         nchroma=0;
         for(int i=0; i<natoms-1; ++i){
            string ss;
            ss.append(to_string(atoms[i].resNum));
            ss.append(atoms[i].resName);
            if(ss.compare(isolabels[j])==0         && 
               atoms[i].atomName.compare("C")==0   &&
               atoms[i+1].atomName.compare("O")==0){
                  chrom_Clist.push_back(i);
                  nchroma++;
               }
         }
         printf("   %s (found %d times)  \n",isolabels[j].c_str(), nchroma);
      }
      nchrom = chrom_Clist.size();
      printf("   Total %d C=O groups are labeled.\n",nchrom);
   }else{
     printf(" Error! Type of spectrum to calculate is not recognized: %s. \n ",jobType.c_str());
     exit(EXIT_FAILURE);
   }

}
