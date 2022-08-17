#include "system.h"

System::System(string gro_file, vector<string> itp_files) 
               : gro_file(gro_file), itp_files(itp_files)
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
   readItp();
   //common();
   //readAtoms();

   // determine starting indices for each residue
   resStart.push_back(0);
   for(unsigned int ii=0; ii<atoms.size()-1; ++ii) 
      if(atoms[ii].resNum != atoms[ii+1].resNum)
         resStart.push_back(ii+1);
}

void System::readItp()
{
// read all itp files provided

   int nItpFiles = itp_files.size();
   printf("\n** Reading GROMACS topology files. %d files provided.**\n",nItpFiles);

   string line, entry, substr, atom_name, residue_name;
   int pos, resid;
   bool recording;

   //int counter=0;
   //int atom_count=0;

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
      resid=-1;
      // 2. Reading the file
      while(getline(itpfile, line)) {
         // Gromacs itp file
         // ignore comment lines starting with ";"
         // and any other line until we find "[ atoms ]"
         // break when we hit an empty line
         //if(line.length()==0) continue;

         pos = line.find(";");
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

         // find [ atoms ] line
         if(line.compare("[ atoms ]")==0){
            recording=true;
            starting=true;
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
            temp.resNum  = stoi(strs[2]);
            temp.resName = strs[3];
            temp.atomName= strs[4];
            temp.chgn    = stoi(strs[5]);
            temp.charge  = stof(strs[6]);
            if(strs.size()==8)
               temp.mass    = stof(strs[7]);

            if(resid==stoi(strs[2]))
            // continue with the same residue
            resid = stoi(strs[2]);

         } // end recording

      } // end reading itp file
      itpfile.close();

   }


}

void System::common()
{
//
// here we detect common repeating molecules like solvent
// currently implemented: water
// 
// 1. check for water; this is good for up to 5-site water
// models but can currently recognize 3 and 4-site models
//   water=false;
//
//   for(int i=0; i<(natoms-4); ++i){
//       if(atoms[i].atomName   =="OW"  &&
//          atoms[i+1].atomName =="HW1" &&
//          atoms[i+2].atomName =="HW2"){
//             water=true;
//             W.sites=3;
//       }
//       if(atoms[i+3].atomName == "MW")
//          W.sites = 4;
//   }
//
//   W.nwater = 0;   
//   if(water){
//      for(int i=0; i<(natoms-W.sites+1); ++i){
//         if(atoms[i].atomName   =="OW"  &&
//            atoms[i+1].atomName =="HW1" &&
//            atoms[i+2].atomName =="HW2"){
//            if(W.sites==3){
//               W.nwater+=1;
//               W.Oind.push_back(i);
//            }
//            else if(W.sites==4){
//               if(atoms[i+3].atomName == "MW"){
//                  W.nwater+=1;
//                  W.Oind.push_back(i);
//               }
//            }
//         }
//      }
//    }
//   
//    if(water)  
//       printf("   %d-site water is detected. Number of water molecules: %d\n",
//              W.sites,W.nwater); 
//
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

void System::readAtoms()
{
//
// Read charges and other info from atoms file or
// itp file (TODO)
//
//   ifstream file(ams_file);
//   if (!file.good()) {
//      printf("ERROR: atom info file %s cannot be read.\n",ams_file.c_str());
//      exit(EXIT_FAILURE);
//   }
//
//   printf("\n** Reading file with atomic information: %s **\n",ams_file.c_str());
//
//   int pos;
//   string line;
//   vector<bool> assigned(natoms,false);
//
//   while(getline(file,line)){
//     
//      pos = line.find("#");
//      if(pos==0) continue;
//
//      string buf;
//      stringstream  linestream(line);
//
//      //linestream >> entry;
//      vector<string> strs;
//      while(linestream >> buf)
//         strs.push_back(buf);
//
//      if(strs.size() != 5){
//         printf("Error! Cannot understand this line in %s file:\n",ams_file.c_str());
//         printf("       %s \n",line.c_str());
//         exit(EXIT_FAILURE);
//      }
//
//      // assign
//      for(int i=0; i<natoms; ++i){
//         if(atoms[i].resNum==stoi(strs[0]) &&
//            atoms[i].resName == strs[1]    && 
//            atoms[i].atomName== strs[2]){
//            atoms[i].charge = stof(strs[3]);
//            atoms[i].mass   = stof(strs[4]);
//            assigned[i] = true;
//         }
//      }
//
//      // we should also check for solvent and other commonly repeating 
//      // residues detected above
//      if(water){
//        // this atom can be a part of water molecule
//        if(strs[2]=="OW")
//            W.Ochg = stof(strs[3]);
//        if(strs[2]=="HW1" || strs[2]=="HW2")
//            W.Hchg = stof(strs[3]);
//        if(strs[2]=="MW")
//            W.Mchg = stof(strs[3]);
//      }
//
//   }
//
//   // if not all atoms have been assigned this might be due
//   // to repeating solvent residues. Check if we have solvent
//   if(water){
//      for(int i=0; i<natoms;++i){
//         if(atoms[i].atomName== "OW"){
//            atoms[i].charge = W.Ochg;
//            atoms[i].mass   = constants::mO;
//            assigned[i] = true;
//         }else if(atoms[i].atomName== "HW1" || atoms[i].atomName=="HW2"){
//            atoms[i].charge = W.Hchg;
//            atoms[i].mass   = constants::mH;
//            assigned[i] = true;
//         }else if(atoms[i].atomName=="MW"){
//            atoms[i].charge = W.Mchg;
//            atoms[i].mass   = 0.0;
//            assigned[i] = true;
//         }
//      }
//   }
//
//   // check again
//   int notas=0;
//   bool all_assigned=true;
//   for(int i=0; i<natoms; ++i){
//      if(!assigned[i]){
//         all_assigned=false;
//         printf("Error! Cannot assign charge and mass to the following atom: %d%s %s\n",
//                 atoms[i].resNum,atoms[i].resName.c_str(),atoms[i].atomName.c_str()); 
//         notas++;
//      }
//   }
//
//   if(!all_assigned){
//      printf("   %d atoms were not assigned charge and mass.\n",notas);
//      exit(EXIT_FAILURE);
//   }
//
//   printf("   All atoms have been assigned charges! Printing them to atoms.info \n");
//
//   ofstream aout_file;
//   aout_file.open("atoms.info");
//   if(!aout_file.is_open()){
//      printf(" Error! Cannot open file: atoms.info \n");
//      exit(EXIT_FAILURE);
//   }
//
//   for(int j=0; j<natoms; ++j)
//      aout_file << atoms[j].resNum << atoms[j].resName << "\t"
//                << atoms[j].atomName << "\t" << atoms[j].charge << "\t" 
//                << atoms[j].mass << endl;  
//
//   aout_file.close();
//
}

string System::exractAndTrim(const string &s, const int a, const int b)
{
   string out = s.substr(a,b);
   remove_leading_trailing_whitespace(out); 
   return out;
}
