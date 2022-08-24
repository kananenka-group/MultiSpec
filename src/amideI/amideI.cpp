#include "amideI.h"

amideI::amideI(string gro_file, string traj_file, vector<string> itp_files,
               string top_file, string spec_type, vector<string> isolabels,
               string nn_map_name, int nframes, int startframe, bool ir,
               float isoShift) :
               s(gro_file, itp_files, top_file), traj_file(traj_file), 
               jobType(spec_type), isolabels(isolabels), map_nn(nn_map_name),
               nframes(nframes), startframe(startframe), ir(ir), 
               isoShift(isoShift)
{
//
// At this point all molecular information has been processed.
// Retrieve it here:
//
   atoms = s.getSystemData();
   natoms= s.getNatoms();   

   // load trajectory file:
   Traj traj(traj_file.c_str());

   // check input data for consistency:
   if (traj.getNatoms() != natoms) {
      printf("ERROR! gro file does not have the same number of atoms as traj file.\n");
      exit(EXIT_FAILURE);
   }

   // locate amide I groups
   findAmideI();

   // find out the type of job requested:
   // either full calculation
   // isotope labeled
   amideIJob();

   // print out what to calculate and set up
   // some flags
   CalcSQuant();

   // allocate memory
   hf.resize(nchrom2,0.0);
   diag_w_nn.resize(nchrom,0.0);
   diag_w_e.resize(nchrom,0.0);
   tdmuf.resize(nchrom*3,0.0);

   // open output files
   houtfile.open("Hamiltonian.bin", ios::binary | ios::out);
   jobfile.open("job.bin", ios::binary | ios::out);
   if(ir) doutfile.open("Dipole.bin", ios::binary | ios::out);

   printf("\n** Generating Excitonic Hamiltonian for %d frames. **\n",nframes);

   // figure out frame to start reading from xtc file
   if(startframe<1){
      printf("   Cannot start from frame %d, setting starting frame to 1.\n",startframe);
      startframe=1;
   }else{
      printf("   Starting from %d frame\n",startframe);
   }

   nframes += (startframe-1);
   int counter=0;

   if(!amdsiso){
   //////////////////////////////////////////////////////////////////// 
   // 
   //  Case 1. Coupled amide I chromophores
   //
   ////////////////////////////////////////////////////////////////////
      while(traj.next()==0 && counter<nframes) {
         counter++;
         if(counter<startframe)
            continue;

         x = traj.getCoords();
         traj.getBox(box);
         nnfs();
          
         updateEx();
      }
   }

}

amideI::~amideI() { }

void amideI::nnfs()
{
//
// Calculate nearest-neighboring frequency shifts
//
   int thisC, thisRes;   

   // remove all angles since this is the first
   // pass for this frame
   //angleID.clear();
   //psi.clear();
   //phi.clear();
   fill_n(diag_w_nn.begin(), nchrom, 0.0);

   for(int i=0; i<nchrom; ++i){
      thisC   = chrom_Clist[i];
      thisRes = atoms[thisC].resNum;
      diag_w_nn[i]  = calc_N_NN_shift(thisC,thisRes,x);
      diag_w_nn[i] += calc_C_NN_shift(thisC,thisRes,x);
      cout << " diagonal shift " << diag_w_nn[i] << endl;
   }
   exit(0);

}

uint amideI::search2(const int st, const int resI, const string &whichtype,
                     const int ResDiff)
{
//
// Search specific atom within a given residue
//
  if (resI>atoms[st].resNum) {  //search forward
    for (int ii=st+1; ii<natoms; ii++)
      if (atoms[ii].atomName.compare(whichtype)==0 && 
          abs(atoms[ii].resNum-atoms[st].resNum)==ResDiff)
        return ii;
  } else {             //search backward
    for (int ii=st-1; ii>=0; ii--)
      if (atoms[ii].atomName.compare(whichtype)==0 && 
          abs(atoms[ii].resNum-atoms[st].resNum)==ResDiff)
        return ii;
  }

  //didn't find type
  return -1;
}


void amideI::calcAngles(const int atomI, const int resI, const rvec *x,
                        float &phi, float &psi)
{
// Steve Strong
// check if already calculated this angle
//  uint ind;
//  for (ind=0; ind<angleID.size(); ind++)
//    if (angleID[ind]==resI)
//      break;
//  if (ind != angleID.size())
//    return ind;

//  angleID.push_back(resI);

  int i1,i2,i3,i4;

  //calculate psi: N-Ca-C-N
  i1=search2(atomI,resI,"N",0);
  i2=search2(atomI,resI,"CA",0);
  i3=atomI;
  i4=search2(atomI,resI+1,"N",2);
  //psi.push_back(calcDihedral(x[i1],x[i2],x[i3],x[i4]));
  cout << " calculating dihedral " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
  psi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);

  // calculate phi: C-N-Ca-C
  i1=search2(atomI,resI-1,"C",1);
  i2=search2(atomI,resI,"N",0);
  i3=search2(atomI,resI,"CA",0);
  i4=atomI;
  cout << " calculating dihedral " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
  phi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);
  //phi.push_back(calcDihedral(x[i1],x[i2],x[i3],x[i4]));

  //return angleID.size()-1;

}

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
   int nchroma=0;
   printf("\n** Reading job type **\n");
   if(jobType=="full"){
      printf("   All peptide groups will be used to calculate spectra.\n");
      isoShift = 0.0;
      //nchrom = namideI;
      chrom_Clist.assign(amideI_Clist.begin(), amideI_Clist.end());
   }else if(jobType=="iso"){
      printf("   Only isotope labeled residues will contribute to the spectra.\n");
      printf("   The following residues are isotope labeled: ");
      for(uint i=0;  i<isolabels.size();++i)
         printf("%s ",isolabels[i].c_str());
      printf("\n");
      // loop over all atoms creating a chromophoe list 
      printf("   Selecting C=O groups: \n");
      for(uint j=0; j<isolabels.size();++j){
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
         printf("   %s (%d)  \n",isolabels[j].c_str(), nchroma);
      }
      printf("   Total %d C=O groups are labeled.\n",nchroma); //chrom_Clist.size());
      printf("   Isotope frequency shift: %7.4f \n",isoShift);
   }else{
     printf(" Error! Type of spectrum to calculate is not recognized: %s. \n ",jobType.c_str());
     exit(EXIT_FAILURE);
   }

   nchrom = (int) chrom_Clist.size();
   nchrom2= nchrom*nchrom;
   printf("   Number of chromphores: %d\n",nchrom);
   //
   // one chrmophore cases are treated separately
   //
   if(nchrom==1) amdsiso=true;
 
   if(nchrom==0){
      printf("Error! Cannot find requested residues.\n");
      exit(EXIT_FAILURE);
   }

}

void amideI::CalcSQuant()
{
   printf("\n** Reading physical properties to be calculated. **\n");
   if(ir){
      printf("   IR flag is on. Transition dipoles will be calculated.\n");
   }else{
      printf("Error! Don't know what to calculate. IR option is available, but not set.\n");
      exit(EXIT_FAILURE);
   }
}

float amideI::calcDihedral(const rvec &x1, const rvec &x2,
                           const rvec &x3, const rvec &x4)
{
  rvec b1,b2,b3;
  addRvec(x2,x1,b1,-1);
  addRvec(x3,x2,b2,-1);
  addRvec(x4,x3,b3,-1);

  pbc(b1,box);
  pbc(b2,box);
  pbc(b3,box);
  
  normalize(b1);
  normalize(b2);
  normalize(b3);
  
  rvec n1,n2;
  cross(b1,b2,n1);
  cross(b2,b3,n2);
  
  //see  https://math.stackexchange.com/a/47084
  rvec m1;
  cross(n1,b2,m1);
  
  float d1,d2;
  d1=dot(n1,n2);
  d2=dot(m1,n2);
  
  //this returns answer in -180 to 180
  float ans = atan2(d2,d1)*180.0/M_PI;
  return ans;
  
  //this returns angle between 0 and 180;
  //return acosd(dot(n1,n2));
}

float amideI::calc_N_NN_shift(const int thisC, const int thisRes, const rvec *x)
{
//
// first check that there is a residue on the N-side
//
   int locC = search2(thisC,thisRes,"C",1);
   if(locC<0)
     return 0.0;

   // if there is then calculate the shift
   float phi=0.;
   float psi=0.;
   //calculate psi: N-Ca-C-N
   int i1,i2,i3,i4;
   i1=search2(thisC,thisRes,"N",0);
   i2=search2(thisC,thisRes,"CA",0);
   i3=thisC; //atomI;
   i4=search2(thisC,thisRes+1,"N",1);
   psi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);
   cout << " N shift atoms psi " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

   // calculate phi: C-N-Ca-C
   i1=search2(thisC,thisRes,"C",1);
   i2=search2(thisC,thisRes,"N",0);
   i3=search2(thisC,thisRes,"CA",0);
   i4=thisC; //atomI;
   phi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);
   cout << " N shift atoms phi " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

   return map_nn.getNNshift(phi,psi,"Nterm");
}


float amideI::calc_C_NN_shift(const int thisC, const int thisRes, const rvec *x)
{
//
// first check that there is a residue on the N-side
//
   int locN = search2(thisC,thisRes+1,"N",2);
   if(locN<0)
     return 0.0;

   // if there is then calculate the shift
   //int nextC = search2(thisC,thisRes+1,"C",1);
   float phi=0.; 
   float psi=0.;
   int i1,i2,i3,i4;
   //calculate psi: N-Ca-C-N
   i1=search2(thisC,thisRes+1,"N",1);
   i2=search2(thisC,thisRes+1,"CA",1);
   i3=search2(thisC,thisRes+1,"C",1);
   i4=search2(thisC,thisRes+1,"N",2);
   psi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);
   cout << " C shift atoms psi " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

   // calculate phi: C-N-Ca-C
   i1=thisC; //search2(thisC,resI-1,"C",1);
   i2=search2(thisC,thisRes+1,"N",1);
   i3=search2(thisC,thisRes+1,"CA",1);
   i4=search2(thisC,thisRes+1,"C",1); //atomI;
   phi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);
   cout << " C shift atoms  phi " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

   //cout << " found next C " << nextC << endl;
   //calcAngles(nextC,thisRes+1,x,phi,psi); 
   //cout <<  " phi and psi " << phi << " " << psi<< endl;
   return map_nn.getNNshift(phi,psi,"Cterm");
}

void amideI::updateEx()
{
// build excitonic Hamiltonian here
   fill_n(hf.begin(), nchrom2, 0.0);
// 1. Diagonal part   
   for(int ii=0; ii<nchrom; ++ii)
      hf[ii*nchrom+ii] = diag_w_nn[ii] + isoShift;

}

