#include "system.h"

System::System(string gro_file, vector<string> itp_files, string top_file) 
               : gro_file(gro_file), itp_files(itp_files), top_file(top_file)
{
// This class will contain information about each atom:
//   To be read from gro file
//   -----------------------------
//   - what residue it belongs to
//   - residue number
//   - atom name
//   To be read from itp (atoms) file
//   ----------------------------
//   - atom charge
//   - atom mass
//

   readGro();
   readTop();
   readItp();

   // determine starting indices for each residue
   resStart.push_back(0);
   for(unsigned int ii=0; ii<atoms.size()-1; ++ii) 
      if(atoms[ii].resNum != atoms[ii+1].resNum)
         resStart.push_back(ii+1);
}

void System::readTop()
{
//
// Here we will read info about molecules from system's
// topology file
//
   ifstream topfile(top_file);
   if(!topfile.good()){
      printf(" ERROR: gromacs file %s cannot be read.\n",top_file.c_str());
      exit(EXIT_FAILURE);
   }
   printf("\n** Reading GROMACS topology file: %s **\n",top_file.c_str());

   string line;
   bool recording=false;
   int pos;

   while(getline(topfile, line)) {

      pos = line.find(";");
      if(pos==0) continue;

      pos = line.find("#");
      if(pos==0) continue;

      if(line.erase(line.find_last_not_of(" \n\r\t")+1).length()==0){
         if(!recording)
            continue;
         else
            break;
      }

      if(line.compare("[ molecules ]")==0){
         recording=true;
         continue;
      }

      if(recording){
         string buf;
         stringstream  linestream(line);
         vector<string> strs;
         while(linestream >> buf)
            strs.push_back(buf);

         molecules.push_back(strs[0]);
         nmol.push_back(stoi(strs[1]));
      }
   }
   topfile.close();  

   printf("   Found the following ``molecules'':\n   ");
   for(unsigned int n=0; n<molecules.size();++n)
      printf("%s(%d) ",molecules[n].c_str(),nmol[n]);
   printf("\n");
}

void System::readItp()
{
// read all itp files provided

   int nItpFiles = itp_files.size();
   printf("\n** Reading GROMACS topology files. %d files provided. **\n",nItpFiles);

   string line, entry, substr, atom_name, mt;
   int pos; 
   int n_this_res;
   bool recording, mtfound, skip;
   vector<AtomsRes> Residue;
   vector<bool> assigned(natoms,false);

   for(int ifile=0; ifile<nItpFiles; ++ifile)
   {
      // 1. Check if this file even exist
      ifstream itpfile(itp_files[ifile]);
      if(!itpfile.good()){
         printf(" ERROR: gromacs file %s cannot be read.\n",itp_files[ifile].c_str());
         exit(EXIT_FAILURE);
      }
      printf("   Reading  %s \n",itp_files[ifile].c_str());   

      recording=false;
      mtfound = false;
      skip=false;
      n_this_res=0;
      //this_residue=-1;
      // 2. Reading the file
      while(getline(itpfile, line)) {
         // Gromacs itp file
         // ignore comment lines starting with ";"
         // and any other line until we find "[ atoms ]"
         // break when we hit an empty line
         //if(line.length()==0) continue;

         pos = line.find(";");
         if(pos==0) continue;

         pos = line.find("#");
         if(pos==0) continue;

         // when searching for [ atoms ] skip empty lines
         // but when recording atoms treat empty lines as section breaks
         // in itp files 
         if(line.erase(line.find_last_not_of(" \n\r\t")+1).length()==0)
         { 
           if(!recording)
              continue;
           else
              break;
         }

         // file [ moleculetype ] first
         if(!mtfound)
            if(line.compare("[ moleculetype ]")==0){
               mtfound = true;
               continue;
            }

         if(mtfound && !skip){
            string buf;
            stringstream  linestream(line);
            vector<string> strs;
            while(linestream >> buf)
               strs.push_back(buf);

            if(strs.size()==2){
               mt = strs[0];
               skip=true;
            }else{
               printf("Error! Cannot understand this moleculetype %s\n",line.c_str());
               exit(EXIT_FAILURE);
            }
         }

         // try to understand what we are dealing with here
         for(unsigned int m=0; m<molecules.size(); ++m)
            if(mt.compare(molecules[m])==0){
               n_this_res = nmol[m];
               continue;
            }

         // find [ atoms ] line
         if(line.compare("[ atoms ]")==0){
            recording=true;
            continue;
         }

         // save atoms 
         if(recording){ 
            string buf;
            stringstream  linestream(line);
            vector<string> strs;
            while(linestream >> buf)
               strs.push_back(buf);
 
            AtomsRes temp;
            temp.resNum  = stoi(strs[2]); //=residue_number
            temp.resName = strs[3];
            temp.atomName= strs[4];
            temp.charge  = stof(strs[6]);

            // sometimes (water) atomic mass is missing
            // trying to add it here
            if(strs.size()>=8){
               temp.mass    = stof(strs[7]);
            }else{
               // probably mass is missing...trying to add
               if(temp.atomName.compare("OW")==0 || 
                  temp.atomName.compare("OH2")==0){
                  temp.mass = 16.00;
               }else if(temp.atomName.compare("HW1")==0 || 
                        temp.atomName.compare("HW2")==0 || 
                        temp.atomName.compare("H1")==0  || 
                        temp.atomName.compare("H2")==0){
                  temp.mass = 1.008;
               }else if(temp.atomName.compare("MW")==0  || 
                        temp.atomName.compare("LP1")==0 || 
                        temp.atomName.compare("LP2")==0){
                  temp.mass = 0.000;
               }else{
                  printf("Error! Not enough atom information in this line: %s \n",line.c_str());
                  exit(EXIT_FAILURE);
               }
            }

            Residue.push_back(temp);

         } // end recording

      } // end reading itp file
      matchRes(Residue,assigned,n_this_res);
      Residue.clear();
      itpfile.close();
   }

   printAtomInfo(assigned);


}

void System::matchRes(vector<AtomsRes> R,vector<bool> &good, int n_this)
{
//
// Match atoms provided in GROMACS configuration file
// with atoms from topology file and update charges and masses
//
   int rlen = R.size();
   bool match=true;

    //
    // if this molecule occurs just once in the config file
    // we will match it by the residue number, residue name, and
    // atom name
    //
    // if this molecule occurs multiple times in the config file
    // we will match it by residue name and atom name
    //
    for(int n=0; n<(natoms-rlen+1); ++n){
       match=true;
       for(int m=0; m<rlen; ++m){
          if(n_this==1){
             if(R[m].resNum!=atoms[n+m].resNum   ||
                R[m].resName!=atoms[n+m].resName || 
                R[m].atomName!=atoms[n+m].atomName){
                match=false;
                continue;
             }
          }else{
             if(R[m].resName!=atoms[n+m].resName ||
                R[m].atomName!=atoms[n+m].atomName){
                match=false;
                continue;
             }
          }
       }
       if(match){
          for(int m=0; m<rlen; ++m){
             atoms[n+m].charge=R[m].charge;
             atoms[n+m].mass=R[m].mass;
             good[n+m]=true;
          } 
       }
    }   

}

void System::readGro()
{
//  From GROMACS manual:
//  This format is fixed, ie. all columns are in a fixed position:
//  - residue number (5 positions, integer)
//  - residue name (5 characters)
//  - atom name (5 characters)
//  - atom number (5 positions, integer)
//
   ifstream file(gro_file);

   if(!file.good()){
      printf(" ERROR: gromacs file %s cannot be read.\n",gro_file.c_str());
      exit(EXIT_FAILURE);
   }

   printf("\n** Reading GROMACS structure file: %s **\n",gro_file.c_str());

   string line, entry, substr;

   int counter=0;
   int atom_count=0;

   AtomsRes temp;

   while(getline(file, line)) {
      counter++;
      stringstream linestream(line);

      // first line is a comment
      if(counter==1)
        continue;     

      // read total number of atoms
      if(counter==2){
         natoms = stoi(line);
         continue;
       }

      // read all atoms:
      substr = exractAndTrim(line,0,5);
      temp.resNum = stoi(substr);

      substr = exractAndTrim(line,5,5);
      temp.resName = substr;

      substr = exractAndTrim(line,10,5);
      temp.atomName = substr;

      atoms.push_back(temp);

      atom_count++;
      counter++;

      if(atom_count==natoms)
        break;
   }

   printf("   Total number of atoms: %d \n",natoms);

}

System::~System(){ }

void System::printAtomInfo(const vector<bool> assigned)
{
   int notas=0;
   bool all_assigned=true;
   for(int i=0; i<natoms; ++i){
      if(!assigned[i]){
         all_assigned=false;
         printf("Error! Cannot assign charge and mass to the following atom: %d%s %s\n",
                 atoms[i].resNum,atoms[i].resName.c_str(),atoms[i].atomName.c_str());
         notas++;
      }
   }

   if(!all_assigned){
     printf("   %d atoms were not assigned charge and mass.\n",notas);
     exit(EXIT_FAILURE);
   }
   
   printf("   All atoms have been assigned charges! Printing them to atoms.info \n");

   ofstream aout_file;
   aout_file.open("atoms.info");
   if(!aout_file.is_open()){
      printf(" Error! Cannot open file: atoms.info \n");
      exit(EXIT_FAILURE);
   }

   for(int j=0; j<natoms; ++j)
      aout_file << atoms[j].resNum << atoms[j].resName << "\t"
                << atoms[j].atomName << "\t" << atoms[j].charge << "\t" 
                << atoms[j].mass << endl;  

   aout_file.close();

}

string System::exractAndTrim(const string &s, const int a, const int b)
{
   string out = s.substr(a,b);
   remove_leading_trailing_whitespace(out); 
   return out;
}
