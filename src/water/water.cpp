#include "water.h"

water::water(string wm_name, string wS_wmap_name, string wB_map_name, 
             string job_type, string traj_file_name, string gro_file_name, 
             string atoms_file_name, int nfr, int d2o,  int nst, bool doir, 
             bool doraman, bool dosfg, bool doDoDov, bool eHp, bool intrac, 
             bool intermc_OHs, float trdip_for_SFG, float fermi_c) : 
             water_model_name_inp(wm_name), jobType(job_type), 
             traj_file(traj_file_name), gro_file(gro_file_name), 
             ams_file(atoms_file_name), wms(wS_wmap_name, job_type, intrac), 
             wmb(wB_map_name, job_type), nframes(nfr), nd2o(d2o), startframe(nst),
             ir(doir), raman(doraman), sfg(dosfg), DoDv(doDoDov), printHam(eHp),
             intermcs(intermc_OHs), tdSFG(trdip_for_SFG), fc(fermi_c)
{
   // removing old files
   vector<string> files_to_remove;
   files_to_remove.push_back("job.bin");
   files_to_remove.push_back("Fz.bin");
   removeFile(files_to_remove);

   // create trajectory object
   Traj traj(traj_file.c_str());

   // read gro file:
   readGro();
   if (traj.getNatoms() != natoms) {
      printf("ERROR! gro file does not have the same number of atoms as traj.\n");
      exit(EXIT_FAILURE);
   }

   // read charges
   readCharges();

   // create water model
   waterModel();

   // and type of job
   waterJob();

   // print out what to calculate and
   // set up some flags there
   CalcSQuant();

   // allocate memory
   nchromt = nchroms + nchromb;
   nchromt2=nchromt*nchromt;

   hf.resize(nchromt2,0.0);
   w10.resize(nchroms,0.0);
   x10.resize(nchroms,0.0);
   p10.resize(nchroms,0.0);
   m10.resize(nchroms,0.0);
   efs.resize(nchroms,0.0);
   efb.resize(nchromb,0.0);
   ROH.resize(nchroms,0.0);
   w20b.resize(nchromb,0.0);
   tdmuf.resize(nchromt*3,0.0);
   plzbf.resize(nchromt*6,0.0);
   fzf.resize(nchromt,0.0);

   // 
   pz = constants::bond_plz_ratio-1.0;

   // open output files... 
   houtfile.open("Hamiltonian.bin", ios::binary | ios::out);
   jobfile.open("job.bin", ios::binary | ios::out);
   if(ir) doutfile.open("Dipole.bin", ios::binary | ios::out);
   if(raman) poutfile.open("Polarizability.bin", ios::binary | ios::out);
   if(sfg) fzoutfile.open("Fz.bin",ios::binary | ios::out);

   vOHa = new rvec[nchroms];
   vOHu = new rvec[nchroms];
   eft  = new rvec[nchroms];

   printf("\n** Generating Excitonic Hamiltonian for %d frames. **\n",nframes);

   // Figure out frame to start reading from xtc file
   if(startframe<1){
      printf("   Cannot start from frame %d, setting starting frame to 1.\n",startframe);
      startframe=1;
   }else{
      printf("   Starting from %d frame\n",startframe);
   }

   nframes += (startframe-1);
   int counter=0;

   if(ws || wf){
   /////////////////////////////////////////////////////////////////////////////
   //                                                                         //
   // Case 1. Stretch fundamental (and bend overtone)                         //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////
      while(traj.next()==0 && counter<nframes) {
         counter++;
         if(counter<startframe)
           continue;

         x = traj.getCoords();
         traj.getBox(box);
         if(moveM) moveMsite();
         calcEf();
         calcWXPM();
         updateEx();
         if(intermcs) intermC();
         if(ir) trDip();
         if(raman) trPol();
         writeH();
         if(ir) writeD(); 
         if(raman) writeP(); 
         if(sfg) writeFz();
      }
   }
   /////////////////////////////////////////////////////////////////////////////
   //                                                                         //
   // Case 2.                                                                 //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////


   // calc. frequency distributions
   if(printHam)   
      freqDist();

   // write a temporary file for passing some info for subsequent job
   //jobfile.write(reinterpret_cast<char*>(&nchromt), sizeof(int));
   writeJ();


}

water::~water() {
  delete [] vOHa;
  delete [] vOHu;
  delete [] eft;
}

void water::CalcSQuant()
{
   printf("\n** Reading physical properties to be calculated. **\n");
   if(ir)
      printf("   IR flag is on. Transition dipoles will be calculated.\n");
   if(raman) 
      printf("   Raman flag is on. Transition polarizabilities will be calculated.\n");
   if(sfg){
      printf("   SFG flag is on. Transition probablities\n");
      printf("                   and scaled transition dipoles will be calculated.\n");
      printf("   See P. A. Pieniazek et al., J. Chem. Phys. 135, 044701 (2011) for more details.\n");                
      printf("   Cut-off in the switching function: %7.5f [A]\n",constants::SWITCHF_CUT);
      if(tdSFG<0){
         printf(" Error: Wrong value of SFG transition dipole location : %7.5f [A] ! \n",tdSFG); 
         exit(EXIT_FAILURE);
      }
      printf("   OH Transition dipole is located %7.5f [A] away from the O atom.\n",tdSFG);
      tdSFG *= constants::A0;
      ir = true;
      raman = true;
   }
}

void water::waterModel()
{
   printf("\n** Reading water model **\n");

   // transform water model name to upper case
   water_model_name = str_toupper(water_model_name_inp);

   // check water model
   if(water_model_name=="TIP4P" || water_model_name=="E3B2"){
      printf("   Water model : TIP4P [ W. L. Jorgensen et al., J. Chem. Phys. 79, 926-935 (1983) ]\n");
      if(atoms_in_mol!=4){
         printf(" Error ! Water model indicated is (TIP4P or E3B2) is inconsistent with *.gro file.\n");
         exit(EXIT_FAILURE);
      }
      moveM = false;
      // set SFG transition dipole
      trdip = 0.67;
      if(tdSFG<0)
         tdSFG = trdip;
      trdip *= constants::A0;
   }else if(water_model_name=="E3B3" || water_model_name=="TIP4P/2005"){
      if (water_model_name=="E3B3")
         printf("   Water model : E3B3 is based on TIP4P/2005 [ Abascal et al. J. Chem. Phys. 123, 234505 (2005) ]\n");
      if (water_model_name=="TIP4P/2005")
         printf("   Water model : TIP4P/2005 [ Abascal et al. J. Chem. Phys. 123, 234505 (2005) ]\n");
      if(atoms_in_mol!=4){
         printf(" Error ! Water model indicated is (TIP4P/2005 or E3B3) is inconsistent with *.gro file.\n");
         exit(EXIT_FAILURE);
      }
      printf("   M-site location will be adjusted to r(OM) = 0.15 A\n"); 
      moveM = true;
      rom = constants::A0*0.15;
      trdip = 0.67;
      if(tdSFG<0)
         tdSFG = trdip;
      trdip *= constants::A0;
   }else if(water_model_name=="SPCE" || water_model_name=="SPC/E" || water_model_name=="SPC"){
      printf("   Water model : SPCE [ Berendsen et al. J. Phys. Chem. 91, 6269 (1987) ]\n");
      if(atoms_in_mol!=3){
         printf(" Error ! Water model indicated is (SPC or SPCE) is inconsistent with *.gro file.\n");
         exit(EXIT_FAILURE);
      }
      moveM = false;
      trdip = 0.58;
      if(tdSFG<0)
         tdSFG = trdip;
      trdip *= constants::A0;
   }else{
      printf(" Error: %s water model is not recognized ! ",water_model_name.c_str()); 
      printf("        Supported water models are: SPC, SPCE, TIP4P, TIP4P/2005, E3B2, E3B3.\n");
      exit(EXIT_FAILURE);
   }
   

}

void water::waterJob(){
/////////////////////////////////////////////////////////////////////////////////////////
//
// Decide what type of job is requested and set up some variables
//
/////////////////////////////////////////////////////////////////////////////////////////
   ws = false;
   wb = false;
   wf = false;
   uncs = false;

   printf("\n** Reading job type ** \n");
   if(jobType=="wsOH"){
      ws = true;
      nchroms = 2*nwater;
      offs = 2;
      printf("   Job Type: pure H2O stretch simulation, %d chromophores.\n",nchroms);
   }else if(jobType=="wsOD"){
      ws = true;
      nchroms = 2*nwater;
      offs = 2;
      printf("   Job Type: pure D2O stretch simulation, %d chromophores.\n",nchroms);
   }else if(jobType=="wsiso"){
      ws = true;
      nchroms = 2*nwater;
      offs = 2;
      printf("   Job Type: Mixed H2O/D2O simulation, %d chromophores.\n",nchroms);
      IsoMix();
   }else if(jobType=="wswbH2O"){
      wf = true;
      nchroms = 2*nwater;
      nchromb = nwater;
      offs = 3;
      printf("   Job Type: pure H2O stretch fundamental-bend overtone simulation.\n"); 
      printf("             %d OH chromophores, %d HOH chromophores \n",nchroms,nchromb);
      printf("             Fermi coupling = %7.2f [cm-1] \n",fc);
   }else if(jobType=="wswbD2O"){
      wf = true;
      nchroms = 2*nwater;
      nchromb = nwater;
      offs = 3;
      printf("   Job Type: pure D2O stretch fundamental-bend overtone simulation.\n");
      printf("             %d OD chromophores, %d DOD chromophores \n",nchroms,nchromb);
      printf("             Fermi coupling = %7.2f [cm-1] \n",fc);
   }else if(jobType=="wswbiso"){
      wf = true;
      nchroms = 2*nwater;
      nchromb = nwater;
      offs = 3;
      printf("   Job Type: Mixed H2O/D2O stretch fundamental-bend overtone simulation.\n");
      printf("             Fermi coupling = %7.2f [cm-1] \n",fc);
      IsoMix();
   }else{
      printf(" Error! Type of spectrum to calculate is not recognized: %s. \n ",jobType.c_str());
      exit(EXIT_FAILURE);
   }

}

void water::readGro(){
///////////////////////////////////////////////////////////////////////////////////////// 
//
// Reading gro file and learning the types of atoms in the system
//
/////////////////////////////////////////////////////////////////////////////////////////
   ifstream file(gro_file);

   if(!file.good()){
      printf(" ERROR: gromacs file %s cannot be read.\n",gro_file.c_str());
      exit(EXIT_FAILURE);
   }

   printf("\n** Reading Gromacs file: %s **\n",gro_file.c_str());

   string line; 
   //bool readFlag=false;

   int counter=0;

   while(getline(file, line)) {

      string entry;
      stringstream  linestream(line);

      if(counter==1){
         linestream >> entry;
         natoms = stoi(entry);
      }else if(counter>1){

         linestream >> entry;

         linestream >> entry;
         aAtoms.push_back(entry);
      }
      counter++;

      if(counter>=(natoms+2))
        continue;
   }

   printf("   Total number of atoms: %d \n",natoms);

   // determine unique atoms
   bool found;
   uAtoms.push_back(aAtoms[0]);
   for(int ii=0; ii<natoms; ++ii){
      int ul = uAtoms.size();
      found = false;
      for(int jj=0; jj<ul; ++jj)
         if(uAtoms[jj]==aAtoms[ii])
           found = true;
      if(!found)
         uAtoms.push_back(aAtoms[ii]);
   }

   printf("   Atoms found: ");
   for(uint ii=0; ii<uAtoms.size(); ++ii)
      printf(" %s ",uAtoms[ii].c_str());
   printf("\n");

   now = 0;
   nhw = 0;
   nmw = 0;

   for(int ii=0; ii<natoms; ++ii){
      if(aAtoms[ii]=="OW")
         now++;
      if(aAtoms[ii]=="HW1" || aAtoms[ii]=="HW2")
         nhw++;
      if(aAtoms[ii]=="MW")
         nmw++;  
   }

   if(now==0){
      printf(" ERROR: No OW atoms are found. This is supposed to be a water simulation!\n");
      exit(EXIT_FAILURE);
   }

   if(nhw==0){
      printf(" ERROR: No HW1 and HW2 atoms are found. This is supposed to be a water simulation!\n");
      exit(EXIT_FAILURE);
   }

   if(nmw>0){
      printf("   MW atoms are found, 4-site water molecule is used.\n"); 
      atoms_in_mol = 4;   
   }else{
      printf("   MW atoms not found, 3-site water molecule is used.\n");
      atoms_in_mol = 3;
   }

   nwater = now;
   printf("   Number of water molecules: %d \n",nwater);   

   if(natoms!=(now+nhw+nmw)){
      printf(" ERROR: The number of OW+HW+MW atoms does not add up \n"); 
      printf("        to the total number of atoms read from %s at \n",gro_file.c_str());
      printf("        this time only pure water can be simulated.\n");
      exit(EXIT_FAILURE);
   }else{
      printf("   It looks like the system is pure water.\n");
   }

   // oxygen atoms
   for(int ii=0; ii<natoms; ++ii)
      if(aAtoms[ii]=="OW")
         oxyInd.push_back(ii);
  
}

void water::readCharges()
{
   ifstream file(ams_file);

   if (!file.good()) {
      printf("ERROR: atom info file %s cannot be read.\n",ams_file.c_str());
      exit(EXIT_FAILURE);
   }

   printf("\n** Reading file with atomic information: %s **\n",ams_file.c_str());

   uChg.resize(uAtoms.size());
   uMas.resize(uAtoms.size());

   string line;

   // read atomic charges and masses line by line from the corresponding file
   int pos;
   while(getline(file,line)){

      pos = line.find("#");
      if(pos==0) continue;

      string buf;
      stringstream  linestream(line);

      //linestream >> entry;
      vector<string> strs;
      while(linestream >> buf)
        strs.push_back(buf);

      // assign charges and masses
      for(uint ii=0; ii<uAtoms.size(); ++ii){
         if(strs[0]==uAtoms[ii]){
            if(strs.size() != 3){
               printf("Error! Wrong line for atom %s in %s. Expecting both charge and mass. \n",uAtoms[ii].c_str(),ams_file.c_str());
               printf(" Expecting : <atom> <charge> <mass>. Got \n");
               for(uint ij=0; ij<strs.size(); ++ij)
                  printf(" %s ",strs[ij].c_str());
               printf("\n");
               exit(EXIT_FAILURE);
            }else{
               uChg[ii] = stof(strs[1]);
               uMas[ii] = stof(strs[2]);
            }
         }
      }
   }

   printf("   Charges found and assigned: ");
   for(uint ii=0; ii<uAtoms.size(); ++ii)
      printf(" %s (%5.3f) ",uAtoms[ii].c_str(),uChg[ii]);
   printf("\n");

   printf("   Masses found and assigned: ");
   for(uint ii=0; ii<uAtoms.size(); ++ii)
      printf(" %s (%5.3f) ",uAtoms[ii].c_str(),uMas[ii]);
   printf("\n");

   aChg.resize(natoms); 
   aMas.resize(natoms);
   
   for(int ii=0; ii<natoms; ++ii)
      for(uint jj=0; jj<uAtoms.size(); ++jj)
         if(aAtoms[ii]==uAtoms[jj]){
            aChg[ii] = uChg[jj];
            aMas[ii] = uMas[jj];
         }

   total_charge = 0.0;
   total_mass   = 0.0;
   for(int ii=0; ii<natoms; ++ii){
      total_charge += aChg[ii];
      total_mass   += aMas[ii];
   }
   inv_total_mass = 1.0/total_mass;

   printf("   Total charge of the system = %5.3f \n",total_charge);
   printf("   Total mass of the system   = %5.3f \n",total_mass);

}

void water::IsoMix()
{
//
   double kst = constants::kEqIsoW;

   printf("\n** Isotope mixed simulation. \n");
   printf("   Equilibrium constant for H2O + D2O <-> 2HOD is %7.5f \n",kst);

   if(nd2o <= 0){
     printf(" Error! Wrong number of D2O molecules: %d \n",nd2o);
     exit(EXIT_FAILURE);
   }

   nh2o = nwater - nd2o;
   printf("   Input: H2O / D2O =  %d / %d \n",nh2o,nd2o);

   double sk, rt, x1, x2;
   int x1i, x2i, xi, ntot;

   sk = sqrt(kst);
   rt = kst*nd2o*nd2o + 16.0*nd2o*nh2o - 2.0*kst*nd2o*nh2o + kst*nh2o*nh2o;

   if(rt >=0 ){
     rt = sqrt(rt);
   }else{
      printf(" Error in solve_isotops ! Sqrt %7.5f \n",rt);
      printf(" Check number of H2O and D2O molecules! \n");
      exit(EXIT_FAILURE);
   }

   x1 = (kst*nd2o + kst*nh2o - sk*rt)/(-8.0 + 2.0*kst);
   x2 = (kst*nd2o + kst*nh2o + sk*rt)/(-8.0 + 2.0*kst);

   x1i = rint(x1);
   x2i = rint(x2);

   if(x1i >= 0 && x1i <= nwater){
      xi = x1i;
   }else if(x2 >= 0 && x2i <= nwater){
      xi = x2i;
   }else{
      printf(" Error! Could not solve an equilibrium problem in IsoMix ! \n");
      exit(EXIT_FAILURE);
   }
   
   // Calculate final number of molecules:
   nh2o = nh2o - xi;
   nd2o = nd2o - xi;
   nhod = 2*xi;
   
   // Make sure we have not lost molecules due to rounding:
   ntot = nh2o + nd2o + nhod;
   if(ntot != nwater){
      printf(" Error in solve_isotops ! Total number of molecules has changed! %d \n",ntot);
      exit(EXIT_FAILURE);
   }

   printf("   Final: H2O / HOD / D2O mixture %d / %d / %d \n",nh2o,nhod,nd2o);

   // determine which molecules are D2O, H2O, and HOD
   vector<int> ovt(oxyInd.size()); 
   iota(ovt.begin(), ovt.end(), 0);

   random_device rd;
   shuffle(ovt.begin(), ovt.end(), rd);

   // assign water molecules to H2O, D2O, and HOD
   mH2O.reserve(nh2o);
   mD2O.reserve(nd2o);
   mHOD.reserve(nhod);

   for(int ii=0; ii<nh2o; ++ii)
      mH2O.push_back(ovt[ii]);

   for(int ii=nh2o; ii<(nh2o+nd2o); ++ii)
      mD2O.push_back(ovt[ii]);
  
   for(int ii=(nh2o+nd2o); ii<(nh2o+nd2o+nhod); ++ii)
      mHOD.push_back(ovt[ii]);

   // for HOD we need to decide which hydroxyl is OH and which is OD
   mt19937 rng(rd());
   uniform_int_distribution<int> uni(1,2);

   for(int ii=0; ii<nhod; ++ii)
      woxyT.push_back(uni(rng));

   // write water indices out
   isofile.open("iso.txt", ios::out);
   isofile << " H2O : ";
   for(int ii=0; ii<nh2o; ++ii)
      isofile << mH2O[ii] << "  ";
   isofile << endl;

   isofile << " D2O : ";
   for(int ii=0; ii<nd2o; ++ii)
      isofile << mD2O[ii] << "  ";
   isofile << endl;;
   
   isofile << " HOD : ";
   for(int ii=0; ii<nhod; ++ii)
      isofile << mHOD[ii] << "  ";
   isofile << endl;
   printf("   H2O / HOD / D2O indices are written to iso.txt file.\n");
}

double water::waterTDC(const rvec &roha, const rvec &trda, const rvec &va, 
                       const rvec &vb, const rvec &box)
{
//
// Transition dipole coupling between 2 hydroxyl chromophores
// of water 
//
   double dm, dm3, wc;
   rvec rohb, trdb, ddv;

   addRvec(va,vb,rohb,-1);
   pbc(rohb,box);
   unitv(rohb);
   addRvec(vb,rohb,trdb,trdip);
   pbc(trdb,box);
   addRvec(trda,trdb,ddv,-1);
   pbc(ddv,box);
   dm = dist(ddv);
   dm3 = dm*dm*dm;
   unitv(ddv);
   wc = constants::AU_TO_WN*(dot(roha,rohb)-3.0*dot(roha,ddv)*dot(rohb,ddv))/dm3;
   return wc;
}

void water::calcWXPM()
{
//
// Diagonal frequency:
//
   int k;
   if(jobType=="wsOH" || jobType=="wsOD" || jobType=="wswbH2O" || jobType=="wswbD2O"){
      for(int i=0;i<nchroms;++i){
         w10[i] = wms.getw01E(efs[i]);
         x10[i] = wms.getx01E(wms.getw01E(efs[i]));
         p10[i] = wms.getp01E(wms.getw01E(efs[i]));
         m10[i] = wms.getm01E(efs[i]);
      }
      for(int k=0; k<nchroms; ++k)
         omegas.push_back(w10[k]);
   }
   else if(jobType=="wsiso" || jobType=="wswbiso"){
      // assign OH
      for(int i=0; i<nh2o; ++i){
         for(int j=0; j<2; ++j){
            k = 2*mH2O[i]+j;
            w10[k] = wms.getw01E_OH(efs[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(efs[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(efs[k]));
            m10[k] = wms.getm01E_OH(efs[k]);
         }
      }
      // assign OD
      for(int i=0; i<nd2o; ++i){
         for(int j=0; j<2; ++j){
            k = 2*mD2O[i]+j;
            w10[k] = wms.getw01E_OD(efs[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(efs[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(efs[k]));
            m10[k] = wms.getm01E_OD(efs[k]);
         }
      }
      // assign HOD
      for(int i=0; i<nhod; ++i){
         if(woxyT[i] == 1){
            k = 2*mHOD[i];
            w10[k] = wms.getw01E_OD(efs[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(efs[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(efs[k]));
            m10[k] = wms.getm01E_OD(efs[k]);
            k = 2*mHOD[i]+1;
            w10[k] = wms.getw01E_OH(efs[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(efs[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(efs[k]));
            m10[k] = wms.getm01E_OH(efs[k]);
         }else{
            k = 2*mHOD[i]+1;
            w10[k] = wms.getw01E_OD(efs[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(efs[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(efs[k]));
            m10[k] = wms.getm01E_OD(efs[k]);
            k = 2*mHOD[i];
            w10[k] = wms.getw01E_OH(efs[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(efs[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(efs[k]));
            m10[k] = wms.getm01E_OH(efs[k]);
         }
      }
      for(int k=0; k<nchroms; ++k)
         omegas.push_back(w10[k]);
   }

   // bend stuff
   if(jobType=="wswbH2O" || jobType=="wswbD2O"){
      for(int i=0;i<nchromb;++i)
         w20b[i] = wmb.getw02E(efb[i]);
   
      for(int k=0; k<nchromb; ++k)
         omegas.push_back(w20b[k]);

   }else if(jobType=="wswbiso"){
      // H2O molecules
      for(int i=0; i<nh2o; ++i){
         k = mH2O[i];
         w20b[k] = wmb.getw01E_HOH(efb[k]) + wmb.getw12E_HOH(efb[k]);
      }
      // D2O molecules
      for(int i=0; i<nd2o; ++i){
         k = mD2O[i];
         if(DoDv){
            w20b[k] = wmb.getw01E_DOD(efb[k]) + wmb.getw12E_DOD(efb[k]);
         }else{
            w20b[k] = 0.0;
         }
      }
      // HOD molecules
      for(int i=0; i<nhod; ++i){
         k = mHOD[i];
         w20b[k] = wmb.getw01E_HOD(efb[k]) + wmb.getw12E_HOD(efb[k]);
      }

      for(int k=0; k<nchromb; ++k)
         omegas.push_back(w20b[k]);
   }


}

void water::intermC(){
///////////////////////////////////////////////////////////////////////
//
//  Calculate intermolecular couplings and update excitonic Hamiltonian
//
///////////////////////////////////////////////////////////////////////

  rvec roha, trda;
  double val;

  for(int i=0; i<nwater; ++i){
     for(int k=1;k<3;++k){
        addRvec(x[oxyInd[i]+k],x[oxyInd[i]],roha,-1);
        pbc(roha,box);
        unitv(roha); 
        addRvec(x[oxyInd[i]],roha,trda,trdip);
        pbc(trda,box); 
        // Hamiltonian
        for(int j=(i+1); j<nwater; ++j){
           for(int l=1;l<3;++l){
              val = waterTDC(roha, trda, x[oxyInd[j]+l], x[oxyInd[j]], box);
              val *= m10[2*i+k-1]*m10[2*j+l-1]*x10[2*i+k-1]*x10[2*j+l-1];
              hf[nchromt*(offs*i+k-1)+offs*j+l-1] = val;
              //hf[nchromt*(offs*i+k-1)+offs*j+l-1] = waterTDC(roha, trda, x[oxyInd[j]+l], x[oxyInd[j]], box);
              //hf[nchromt*(offs*i+k-1)+offs*j+l-1] *= m10[2*i+k-1]*m10[2*j+l-1]*x10[2*i+k-1]*x10[2*j+l-1];
              w_inter.push_back(val);
           }
        }
     }
  }
  
}


void water::trDip()
{
   rvec out, vcom, slabv, loctd;
   float scl;
   fill_n(tdmuf.begin(), 3*nchromt, 0.0);

   if(sfg) com(vcom);

   for(int n=0; n<nwater; ++n){
      for(int m=0; m<2; ++m){
         scl = m10[2*n+m]*x10[2*n+m];
         sclRvec(vOHa[2*n+m], out, scl);
         if(sfg){
            addRvec(x[oxyInd[n]],vOHa[2*n+m],loctd,tdSFG);
            addRvec(loctd,vcom,slabv,-1.0);
            pbc(slabv,box);
            fzf[offs*n+m] = switchf(ecRvec(slabv,2));
         }
         copyRRvec(&tdmuf[3*(offs*n+m)], out);
      }
   }

   // Zero strength for bend overtones 
   if(wf)
      for(int n=0; n<nwater; ++n)
         tdmuf[3*(offs*n+2)] = 0.0;
   
}

void water::trPol()
{
   vector<float> out(6,0.0);
   fill_n(plzbf.begin(), 6*nchromt, 0.0);

   for(int n=0; n<nwater; ++n){
      for(int m=0; m<2; ++m){
         OuterRvec(&out[0], vOHa[2*n+m]);
         plzbf[6*(offs*n+m)]   = (pz*out[0] + 1.0)*x10[2*n+m];
         plzbf[6*(offs*n+m)+1] = pz*out[1]*x10[2*n+m];
         plzbf[6*(offs*n+m)+2] = pz*out[2]*x10[2*n+m];
         plzbf[6*(offs*n+m)+3] = (pz*out[3] + 1.0)*x10[2*n+m];
         plzbf[6*(offs*n+m)+4] = pz*out[4]*x10[2*n+m];
         plzbf[6*(offs*n+m)+5] = (pz*out[5] + 1.0)*x10[2*n+m];
      }
   }

   // Zero strength for bend overtone
   if(wf){
      for(int n=0; n<nwater; ++n){
         plzbf[6*(offs*n+2)]   = 0.0;
         plzbf[6*(offs*n+2)+1] = 0.0;
         plzbf[6*(offs*n+2)+2] = 0.0;
         plzbf[6*(offs*n+2)+3] = 0.0;
         plzbf[6*(offs*n+2)+4] = 0.0;
         plzbf[6*(offs*n+2)+5] = 0.0;
      }
   }

}

void water::writeH()
{
  int jj = 0;
  int pr_size;
  pr_size = nchromt;
  for(int ii=0; ii<nchromt; ++ii){
     jj = ii*nchromt+ii;
     houtfile.write(reinterpret_cast<char*>(&hf[jj]), pr_size*sizeof(float));
     pr_size -= 1;
  }
}

void water::writeJ()
{ jobfile.write(reinterpret_cast<char*>(&nchromt), sizeof(int)); } 

void water::writeD()
{ doutfile.write(reinterpret_cast<char*>(&tdmuf[0]), (3*nchromt)*sizeof(float)); }

void water::writeP()
{ poutfile.write(reinterpret_cast<char*>(&plzbf[0]), (6*nchromt)*sizeof(float)); }

void water::writeFz()
{ fzoutfile.write(reinterpret_cast<char*>(&fzf[0]), (nchromt)*sizeof(float)); }

void water::calcEf(){
/////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate electric field on all H atoms; for OH stretch and HOH bend
//
/////////////////////////////////////////////////////////////////////////////////////////
   int eHl, tagH, Ok, thisA;
   rvec roh, rohs, rAH;
   rvec vp, tmp1, tmp2;
   float doh, dah, dah3, scla;

   fill_n(ROH.begin(), nchroms, 0.0);
   fill_n(efs.begin(), nchroms, 0.0);

   for(int n=0; n<nchroms; ++n)
      setRvec(eft[n], 0.0);

   for(int i=0;i<nwater;++i){
      for(int j=1;j<3;++j){
         eHl = 2*i+j-1;
         tagH=oxyInd[i]+j;
         addRvec(x[tagH],x[oxyInd[i]],rohs,-1);
         pbc(rohs,box);
         ROH[eHl] = dist(rohs);
         copyRvec(rohs,vOHu[eHl]);
         unitv(rohs);
         copyRvec(rohs,vOHa[eHl]);
         for(int k=0;k<nwater;++k){
            if(i==k) continue;
            Ok = oxyInd[k];
            addRvec(x[Ok],x[tagH],roh,-1);
            pbc(roh,box);
            doh = dist(roh);
            if(doh>constants::O_to_H_dist_water_cutoff) continue;
            for(int l=0; l<atoms_in_mol; ++l){
               thisA = oxyInd[k]+l;
               addRvec(x[tagH],x[thisA],rAH,-1);
               pbc(rAH,box);
               dah = dist(rAH);
               dah3 = dah*dah*dah;
               scla = aChg[thisA]/dah3;
               addRvec(rAH, eft[eHl], scla);
            }
         }
         efs[eHl] = dot(eft[eHl], rohs);
      }
   }
   
   if(wb || wf){
      fill_n(efb.begin(), nchromb, 0.0);

      for(int n=0; n<nwater; ++n){
         cross(vOHu[2*n], vOHu[2*n+1], vp);
         unitv(vp);      
         cross(vp, vOHu[2*n], tmp1);
         cross(vOHu[2*n+1], vp, tmp2);
         efb[n] = dot(tmp1, eft[2*n])/ROH[2*n] + dot(tmp2, eft[2*n+1])/ROH[2*n+1]; 
      }
   }

}

void water::updateEx(){
//////////////////////////////////////////////////////////////////////////////////////
//
// Add diagonal elements of exciton Hamiltonian and intramolecular couplings 
//
//////////////////////////////////////////////////////////////////////////////////////

  int ij;
  double val;

  fill_n(hf.begin(), nchromt2, 0.0);

  // diagonal...
  for(int ii=0; ii<nwater; ++ii){
     ij = offs*ii;
     hf[ij*nchromt+ij] = w10[2*ii];
     ij = offs*ii + 1;
     hf[ij*nchromt+ij] = w10[2*ii+1];
  }

  // bend overtone
  if(wf)
     for(int ii=0; ii<nwater; ++ii)
        hf[nchromt*(3*ii+2) + 3*ii+2]  = w20b[ii];

  // intramolecular coupling...
  for(int ii=0; ii<nwater; ++ii){
     ij = offs*ii;
     val = wms.getcnn(efs[2*ii],x10[2*ii],p10[2*ii], efs[2*ii+1],x10[2*ii+1],p10[2*ii+1]);
     hf[ij*nchromt+ij+1] = val;

     //hf[ij*nchromt+ij+1] = wms.getcnn(efs[2*ii],x10[2*ii],p10[2*ii],
     //                                 efs[2*ii+1],x10[2*ii+1],p10[2*ii+1]);
     w_intra.push_back(val);
  }

  // Fermi coupling
  if(wf){ 
     for(int ii=0; ii<nwater; ++ii){
        hf[3*ii*nchromt+3*ii+2] = fc;
        hf[(3*ii+1)*nchromt+(3*ii+1)+1] = fc;
     }
  }

}

void water::com(rvec &vcom)
{
// determine center-of-mass of a slab.
   rvec out;
   setRvec(out, 0.0);

   for(int ii=0; ii<nwater; ++ii)
      for(int l=0; l<atoms_in_mol; ++l)
         addRvec(x[oxyInd[ii]+l], out, aMas[atoms_in_mol*ii+l]);

   sclRvec(out, vcom, inv_total_mass);
}

void water::moveMsite()
{
   int tagM;
   rvec om;

   for(int ii=0; ii<nwater; ++ii){
      tagM=oxyInd[ii]+3;
      addRvec(x[tagM],x[oxyInd[ii]],om,-1);
      pbc(om,box);
      unitv(om);
      addRvec(x[oxyInd[ii]],om,x[tagM],rom);
   }
}

void water::freqDist()
{
   uint t;
//
// Diagonal frequencies
//
   if(!DoDv)
      omegas.erase(remove(begin(omegas), end(omegas), 0.0), end(omegas));

   ofstream w_dist_file;
   w_dist_file.open(w_dist_fname);
   if(!w_dist_file.is_open()){
      printf(" Error! Cannot open file: %s \n",w_dist_fname.c_str());
      exit(EXIT_FAILURE);
   }

   printf("   Writing diagonal frequencies into %s \n",w_dist_fname.c_str());
   for(t=0; t<omegas.size(); ++t)
      w_dist_file << omegas[t] << endl;

   w_dist_file.close();
//
// Intramolecular couplings
//
   ofstream w_intra_file;
   w_intra_file.open(w_intra_fname);
   if(!w_intra_file.is_open()){
      printf(" Error! Cannot open file: %s \n",w_intra_fname.c_str());
      exit(EXIT_FAILURE);
   }

   printf("   Writing intramolecular couplings [in cm-1] into %s \n",w_intra_fname.c_str());
   if(jobType=="wswbH2O" || jobType=="wswbD2O" || jobType=="wswbiso")
      printf("   Note that Fermi couplings %7.2f [cm-1] will not be printed.\n",fc);

   for(t=0; t<w_intra.size(); ++t)
      w_intra_file << w_intra[t] << endl;

   w_intra_file.close();
//
// Intermolecular couplings
//
  ofstream w_inter_file;
   w_inter_file.open(w_inter_fname);
   if(!w_inter_file.is_open()){
      printf(" Error! Cannot open file: %s \n",w_inter_fname.c_str());
      exit(EXIT_FAILURE);
   }

   printf("   Writing intermolecular couplings [in cm-1] into %s \n",w_inter_fname.c_str());
   for(t=0; t<w_inter.size(); ++t)
      w_inter_file << w_inter[t] << endl;

   w_inter_file.close();
   
}
