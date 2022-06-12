#include "water.h"

water::water(string wm_name, string ctype, int nfr, string wS_wmap_name, string job_type,
             string traj_file_name, string gro_file_name, string charge_file_name, 
             bool doir, bool doraman, bool dosfg, int d2o) :
             water_model_name(wm_name), chromType(ctype), nframes(nfr), 
             jobType(job_type), traj_file(traj_file_name), gro_file(gro_file_name), 
             chg_file(charge_file_name), wms(wS_wmap_name, job_type),
             ir(doir), raman(doraman), sfg(dosfg), nd2o(d2o)
{
   // create trajectory object
   Traj traj(traj_file.c_str());

   // read gro file:
   readGro();

   // read charges
   readCharges();

   // create water model
   waterModel();

   // find out the type of chromphore
   waterChrom(); 

   // and type of job
   waterJob();

   rvec roha, trda; 

   int ii, pp;

   ef.assign(nchrom,0.0);

   ndim = nchrom*(nchrom+1)/2;
   hf.resize(ndim);
   w10.resize(nchrom);
   x10.resize(nchrom);
   p10.resize(nchrom);
   m10.resize(nchrom);

   // transition dipoles
   tdmuf.resize(nchrom*3);

   // transition polarizabilities
   plzbf.resize(nchrom*6);
   pz = bond_plz_ratio-1.0;

   printf("\n** Generating Excitonic Hamiltonian for %d frames. **\n",nframes);
   int counter=0;
   // 
   // case 1. Water hydroxyl stretch
   //
   if(ws){
      while(traj.next()==0 && counter<nframes) {
         x = traj.getCoords();
         traj.getBox(box);
         calcEf();
         calcWXPM();
         updateEx();
         // TD and TDC coupling
         ii=0;
         for(int i=0; i<nwater; ++i){
            for(int k=1;k<3;++k){
               ii += 3-k;
               addRvec(x[oxyInd[i]+k],x[oxyInd[i]],roha,-1);
               pbc(roha,box);
               unitv(roha); 
               addRvec(x[oxyInd[i]],roha,trda,trdip);
               pbc(trda,box); 
               // transition dipole
               if(ir){
                  for(int kk=0; kk<3;++kk)
                     tdmuf[3*(2*i+k-1)+kk] = roha[kk]*m10[2*i+k-1]*x10[2*i+k-1];
               }
               // transition polarizability
               if(raman){
                  pp = 0;
                  for(int kk=0; kk<3;++kk){
                     plzbf[6*(2*i+k-1)+pp] = (pz*roha[kk]*roha[kk]+1.0)*x10[2*i+k-1]; 
                     pp+=1;
                     for(int ll=(kk+1); ll<3; ++ll){
                        plzbf[6*(2*i+k-1)+pp] = pz*roha[kk]*roha[ll]*x10[2*i+k-1]; 
                        pp+=1;
                     }
                  }
               }
               // Hamiltonian
               for(int j=(i+1); j<nwater; ++j){
                  for(int l=1;l<3;++l){
                     hf[ii] = waterTDC(roha, trda, x[oxyInd[j]+l], x[oxyInd[j]], box);
                     hf[ii] *= m10[2*i+k-1]*m10[2*j+l-1]*x10[2*i+k-1]*x10[2*j+l-1];
                     ii++;
                  }
               }
            }
         }
  
         // update excitonic Hamiltonian and dipole derivtv
         ht.insert(ht.end(), begin(hf), end(hf));
         if(ir) tdmut.insert(tdmut.end(), begin(tdmuf), end(tdmuf));
         if(raman) plzbt.insert(plzbt.end(), begin(plzbf), end(plzbf));
         counter++;
      }
   }

   printf("\n** Writing Excitonic Hamiltonian into Hamiltonian.bin **\n");
   ofstream houtfile;
   houtfile.open("Hamiltonian.bin", ios::binary | ios::out);
   houtfile.write(reinterpret_cast<char*>(&ht[0]), ht.size()*sizeof(float));
   houtfile.close();


   if(ir){
      printf("\n** Writing transition dipole moment trajectory into Dipole.bin **\n");
      ofstream doutfile;
      doutfile.open("Dipole.bin", ios::binary | ios::out);
      doutfile.write(reinterpret_cast<char*>(&tdmut[0]), tdmut.size()*sizeof(float));
      doutfile.close();
   }

   if(raman){
      printf("\n** Writing transition polarizability trajectory into Polarizability.bin **\n");
      ofstream poutfile;
      poutfile.open("Polarizability.bin", ios::binary | ios::out);
      poutfile.write(reinterpret_cast<char*>(&plzbt[0]), plzbt.size()*sizeof(float));
      poutfile.close();
   }

}

void water::waterModel()
{
   printf("\n** Reading water model **\n");

   // check water model
   if(water_model_name=="TIP4P" || water_model_name=="tip4p"){
      printf("   Water model : TIP4P [ W. L. Jorgensen et al., J. Chem. Phys. 79, 926-935 (1983) ]\n");
      water_model_name_caps = "TIP4P";
      atoms_in_mol_wm = 4;
      trdip = 0.67*A0;
   }else{
      // SPC/E trdip = 0.58 Angstrom
      printf(" Error: %s water model is not recognized ! ",water_model_name.c_str()); 
      exit(EXIT_FAILURE);
   }
   

}

void water::waterChrom()
{
   ws = false;
   wb = false;
   wf = false;

   printf("\n** Reading chromophore type **\n");

   if(chromType=="ws"){
      printf("   Chromophore: water hydroxyl stretch \n");
      ws = true;
      nchrom = 2*nwater;
   }else if(chromType=="wb"){
      printf("   Chromophore: water bend \n");
      wb = true;
      nchrom = nwater;
   }else if(chromType=="ws+wb"){
      printf("   Chromophore: water hydroxyl stretch + bend (with Fermi resonance)\n");
      wf = true;
      nchrom = 3*nwater;
   }else{
      printf(" Error: %s type of chromophore is not recognized ! \n",chromType.c_str());
      exit(EXIT_FAILURE);
   }
}

void water::waterJob()
{
   printf("\n** Reading job type ** \n");
   if(jobType=="wsOH"){
      pure = true;
      printf("   Job Type: pure H2O stretch simulation, %d chromophores.\n",nchrom);
   }else if(jobType=="wsOD"){
      pure = true;
      printf("   Job Type: pure D2O stretch simulation, %d chromophores.\n",nchrom);
   }else if(jobType=="wsiso"){
      iso = true;
      printf("   Job Type: Mixed H2O/D2O simulation, %d chromophores.\n",nchrom);
      IsoMix();
   }else{
      printf(" Error! Type of spectrum to calculate is not recognized: %s. \n ",jobType.c_str());
      exit(EXIT_FAILURE);
   }


}

void water::readGro()
{
/* 
 *
 * Reading gro file and learning the types of atoms in the system
 *
 */
   ifstream file(gro_file);

   if (!file.good()) {
      printf(" ERROR: gromacs file %s cannot be read.\n",gro_file.c_str());
      exit(EXIT_FAILURE);
   }

   printf("\n** Reading Gromacs file: %s **\n",gro_file.c_str());

   string line; 
   bool readFlag=false;

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
   for(unsigned int ii=0; ii<uAtoms.size(); ++ii)
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
      printf(" ERROR: The number of OW+HW+MW atoms does not add up \n to the total number of atoms read from %s \n at this time only pure water can be simulated.\n",gro_file.c_str());
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
   ifstream file(chg_file);

   if (!file.good()) {
      printf("ERROR: charge file %s cannot be read.\n",chg_file.c_str());
      exit(EXIT_FAILURE);
   }

   printf("\n** Reading charge file: %s **\n",chg_file.c_str());

   uChg.resize(uAtoms.size());

   string line;

   while(getline(file, line)) {

      string entry;
      stringstream  linestream(line);

      linestream >> entry;

      for(unsigned int ii=0; ii<uAtoms.size(); ++ii){
         if(entry==uAtoms[ii]){
            linestream >> entry;
            uChg[ii] = stof(entry);
         }
      }
   }

   printf("   Charges found and assigned: ");
   for(unsigned int ii=0; ii<uAtoms.size(); ++ii)
      printf(" %s (%5.3f) ",uAtoms[ii].c_str(),uChg[ii]);
   printf("\n");

   aChg.resize(natoms); 
   
   for(int ii=0; ii<natoms; ++ii)
      for(unsigned int jj=0; jj<uAtoms.size(); ++jj)
         if(aAtoms[ii]==uAtoms[jj])
            aChg[ii] = uChg[jj];

   total_charge = 0.0;
   for(int ii=0; ii<natoms; ++ii)
      total_charge += aChg[ii];

   printf("   Total charge of the system = %5.3f \n",total_charge);

}

void water::updateEx()
{
  hf.assign(ndim,0.0);   

  // diagonal
  int jj=0;
  for(int ii=0; ii<nchrom; ++ii){
     hf[jj] = w10[ii];
     jj += nchrom-ii;
  }

  // NN coupling
  jj = 1; 
  for(int ii=0; ii<nchrom; ii+=2){
     hf[jj] = wms.getcnn(ef[ii],x10[ii],p10[ii],ef[ii+1],x10[ii+1],p10[ii+1]);
     jj += 2*(nchrom-ii)-1;
  }
  
}

void water::IsoMix()
{
//
   printf("   Isotope mixed simulation. \n");
   printf("   Equilibrium constant for H2O + D2O <-> 2HOD is %7.5f \n",kEqIsoW);

   if(nd2o <= 0){
     printf(" Error! Wrong number of D2O molecules: %d \n",nd2o);
     exit(EXIT_FAILURE);
   }

   nh2o = nwater - nd2o;
   printf("   Input: H2O / D2O =  %d / %d \n",nh2o,nd2o);

   double sk, rt, x1, x2;
   int x1i, x2i, xi, ntot;

   sk = sqrt(kEqIsoW);
   rt = kEqIsoW*nd2o*nd2o + 16.0*nd2o*nh2o - 2.0*kEqIsoW*nd2o*nh2o + kEqIsoW*nh2o*nh2o;

   if(rt >=0 ){
     rt = sqrt(rt);
   }else{
      printf(" Error in solve_isotops ! Sqrt %7.5f \n",rt);
      printf(" Check number of H2O and D2O molecules! \n");
      exit(EXIT_FAILURE);
   }

   x1 = (kEqIsoW*nd2o + kEqIsoW*nh2o - sk*rt)/(-8.0 + 2.0*kEqIsoW);
   x2 = (kEqIsoW*nd2o + kEqIsoW*nh2o + sk*rt)/(-8.0 + 2.0*kEqIsoW);

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
   wc = au_to_wn*(dot(roha,rohb)-3.0*dot(roha,ddv)*dot(rohb,ddv))/dm3;
   return wc;
}

void water::calcEf()
{
//
// calculate electric field on all H atoms of water molecules
//
   int eHl, tagH, Ok, thisA;
   rvec roh, rohs, rAH;
   float doh, dah, dah3;

   vector<float> eft;
   eft.assign(3,0.0);

   for(int i=0;i<nwater;++i){
       for(int j=1;j<3;++j){
          eHl = 2*i+j-1;
          tagH=oxyInd[i]+j;
          addRvec(x[tagH],x[oxyInd[i]],rohs,-1);
          pbc(rohs,box);
          unitv(rohs);
          fill(eft.begin(), eft.end(), 0.0);
          for(int k=0;k<nwater;++k){
             if(i==k) continue;
             Ok = oxyInd[k];
             addRvec(x[Ok],x[tagH],roh,-1);
             pbc(roh,box);
             doh = dist(roh);
             if(doh>O_to_H_dist_water_cutoff) continue;
             for(int l=0; l<atoms_in_mol; ++l){
                thisA = oxyInd[k]+l;
                addRvec(x[tagH],x[thisA],rAH,-1);
                pbc(rAH,box);
                dah = dist(rAH);
                dah3 = dah*dah*dah;
                for(int m=0; m<3; ++m)
                   eft[m] += rAH[m]*aChg[thisA]/dah3;
             }
          }
          ef[eHl] = 0.0;
          for(int n=0; n<3; ++n)
             ef[eHl] += eft[n]*rohs[n];
       }
   }
}

void water::calcWXPM()
{
   int k;
   if((jobType=="wsOH") || (jobType=="wsOD")){
      for(int i=0;i<nchrom;++i){
         w10[i] = wms.getw01E(ef[i]);
         x10[i] = wms.getx01E(wms.getw01E(ef[i]));
         p10[i] = wms.getp01E(wms.getw01E(ef[i]));
         m10[i] = wms.getm01E(ef[i]);
      }
   }
   else if(jobType=="wsiso"){
      // assign OH
      for(int i=0; i<nh2o; ++i){
         for(int j=0; j<2; ++j){
            k = 2*mH2O[i]+j;
            w10[k] = wms.getw01E_OH(ef[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(ef[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(ef[k]));
            m10[k] = wms.getm01E_OH(ef[k]);
         }
      }
      // assign OD
      for(int i=0; i<nd2o; ++i){
         for(int j=0; j<2; ++j){
            k = 2*mD2O[i]+j;
            w10[k] = wms.getw01E_OD(ef[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(ef[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(ef[k]));
            m10[k] = wms.getm01E_OD(ef[k]);
         }
      }
      // assign HOD
      for(int i=0; i<nhod; ++i){
         if(woxyT[i] == 1){
            k = 2*mHOD[i];
            w10[k] = wms.getw01E_OD(ef[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(ef[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(ef[k]));
            m10[k] = wms.getm01E_OD(ef[k]);
            k = 2*mHOD[i]+1;
            w10[k] = wms.getw01E_OH(ef[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(ef[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(ef[k]));
            m10[k] = wms.getm01E_OH(ef[k]);
         }else{
            k = 2*mHOD[i]+1;
            w10[k] = wms.getw01E_OD(ef[k]);
            x10[k] = wms.getx01E_OD(wms.getw01E_OD(ef[k]));
            p10[k] = wms.getp01E_OD(wms.getw01E_OD(ef[k]));
            m10[k] = wms.getm01E_OD(ef[k]);
            k = 2*mHOD[i];
            w10[k] = wms.getw01E_OH(ef[k]);
            x10[k] = wms.getx01E_OH(wms.getw01E_OH(ef[k]));
            p10[k] = wms.getp01E_OH(wms.getw01E_OH(ef[k]));
            m10[k] = wms.getm01E_OH(ef[k]);
         }
      }
   }

}
