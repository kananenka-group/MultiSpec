#include "amideI.h"

amideI::amideI(string gro_file, string traj_file, string mapsFile, 
               vector<string> itp_files,
               string top_file, string spec_type, vector<string> isolabels,
               string nn_map_name, string nnc_map_name, string el_map_name, 
               int nframes, int startframe, bool ir, float isoShift) :
               s(gro_file, itp_files, top_file), traj_file(traj_file), 
               jobType(spec_type), isolabels(isolabels),
               map_nn(nn_map_name), map_nnc(nnc_map_name), 
               map_el(el_map_name,mapsFile), nframes(nframes), 
               startframe(startframe), ir(ir), isoShift(isoShift)
{
//
// At this point all molecular information has been processed.
// Retrieve it here:
//
   atoms  = s.getSystemData();
   natoms = s.getNatoms();   
   chgSt  = s.getChgSt();
   grpInd = s.getCgInd();

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
   tdmuf.resize(nchrom*3,0.0);
   diag_w_nn.resize(nchrom,0.0);
   diag_w_el.resize(nchrom,0.0);

   vecCO = new rvec[nchrom];
   vecCN = new rvec[nchrom];
   trdipL= new rvec[nchrom];
   dip   = new rvec[nchrom];

   ncog  = (int) chgSt.size()-1; 
   cog   = new rvec[ncog];

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
         if(ir) trDip();
         getCCG();
         elstCG();
         //elst();
         updateExd();
         couplings();
         writeH();
         if(ir) writeD();
      }
   }

   // write out files.
   writeJ();

}

void amideI::getCCG()
{
//
// calculate geometric centers of each charge group
//
   for(int ii=0; ii<ncog; ++ii)
      setRvec(cog[ii],0.0);  

   int jj,nthis,thisSt;
   rvec tmpcog,tmpdiff,tmpvec,xref;

   for (int ii=0; ii<ncog; ++ii) {
      thisSt=chgSt[ii];
      nthis=chgSt[ii+1]-thisSt;
      copyRvec(x[thisSt],tmpcog);
      copyRvec(tmpcog,xref);
      for(jj=1; jj<nthis; ++jj){
         copyRvec(x[thisSt+jj],tmpvec);
         addRvec(tmpvec,xref,tmpdiff,-1.0);
         pbcOther(tmpvec,tmpdiff,box);   //get image closest to 1st point
         addRvec(tmpvec,tmpcog,1.0);
      }
      multRvec(tmpcog,1.0/nthis);
      copyRvec(tmpcog,cog[ii]);
   }
}

void amideI::couplings()
{
//
// Calculate NN couplings and TDC here
//
   rvec tmpvec;
   float id, id3, tmpcoup;
   int iC, jC, res1, res2, resdiff;

   for(int i=0; i<nchrom; ++i){
      for(int j=(i+1); j<nchrom; ++j){
         iC = chrom_Clist[i];
         jC = chrom_Clist[j];
         res1 = atoms[iC].resNum;
         res2 = atoms[jC].resNum;
         resdiff = res1 - res2;
         // NN coupling
         // here implement more rigorous check for example
         // two res in ion channels diff S1-S4 will have resdiff+-1
         // but they won't be nearest neighbors
         if(resdiff==1){
            tmpcoup = calc_N_NN_CS(iC,res1,x,2);
            //cout << " nn coup " << res1 << " " << res2 << endl;
         }else if(resdiff==-1){
            tmpcoup = calc_C_NN_CS(iC,res1,x,2);
            //cout << " nn coup " << res1 << " " << res2 << endl;
         }else{
            // TDC coupling
            addRvec(trdipL[i],trdipL[j],tmpvec,-1);
            pbc(tmpvec,box);
            id=1.0/dist(tmpvec);
            id3=id*id*id;

            tmpcoup  = dot(dip[i],dip[j])*id3;
            tmpcoup -= 3.0*dot(dip[i],tmpvec)*dot(dip[j],tmpvec)*id3*id*id;
            tmpcoup *= constants::amideIcoupConst;
         }
         hf[i*nchrom+j] = tmpcoup;
         hf[j*nchrom+i] = tmpcoup;
      }
   }

}

vector<int> amideI::get_excludes(int cid, int type)
{
   vector<int> excl;
   int thisO, thisN, thisCa, thisH;
   int nextO, nextC, nextN, nextCa, nextH;
   int prevO, prevC, prevN, prevCa, prevH;
   int tmpex;

   // Exclude this peptide group and Ca atom of the residue
   // C
   if(type==0){
      excl.push_back(cid);
   }else{
      tmpex = grpInd[cid];
      if(find(excl.begin(), excl.end(), tmpex)==excl.end())
         excl.push_back(tmpex);
   }
   // O
   thisO = search2(cid, atoms[cid].resNum+1, "O", 0);
   if(thisO>0){
      if(type==0){
         excl.push_back(thisO);
      }else{
         tmpex = grpInd[thisO];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[thisO]);
      }
   }
   // N 
   thisN = search2(cid, atoms[cid].resNum+1, "N", 1);
   if(thisN>0){
      if(type==0){
         excl.push_back(thisN);
      }else{
         tmpex = grpInd[thisN];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[thisN]); 
      }
   }
   // H
   thisH = search2(cid, atoms[cid].resNum+1, "H", 1);
   if(thisH>0){
      if(type==0){
         excl.push_back(thisH);
      }else{
         tmpex = grpInd[thisH];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[thisH]);
      }
   }
   // Ca
   thisCa = search2(cid, atoms[cid].resNum+1, "CA", 1);
   if(thisCa>0){
      if(type==0){
         excl.push_back(thisCa);
      }else{
         tmpex = grpInd[thisCa];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[thisCa]);
      }
   }
   // next-neighbor
   // exclude NN peptide group
   nextO = search2(cid, atoms[cid].resNum+1, "O", 1);
   if(nextO>0){
     if(type==0){
        excl.push_back(nextO);
     }else{
        tmpex = grpInd[nextO];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[nextO]);
     } 
   }
   // check if all these atoms are found 
   // at the same time only that way we can be sure it
   // is a backbone group of compare "C" index with 
   // "C" atoms in amideI_list
   nextC = search2(cid, atoms[cid].resNum+1, "C", 1);
   if(nextC>0){
     if(type==0){
        excl.push_back(nextC);
     }else{
        tmpex = grpInd[nextC];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[nextC]);
     }
   }

   nextN = search2(cid, atoms[cid].resNum+1, "N", 2);
   if(nextN>0){
     if(type==0){
        excl.push_back(nextN);
     }else{
        tmpex = grpInd[nextN];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[nextN]);
     }
   }

   nextCa = search2(cid, atoms[cid].resNum+1, "CA", 2);
   if(nextCa>0){
     if(type==0){
        excl.push_back(nextCa);
     }else{
        tmpex = grpInd[nextCa];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[nextCa]);
     }
   }

   nextH = search2(cid, atoms[cid].resNum+1, "H", 2);
   if(nextH>0){
     if(type==0){
        excl.push_back(nextH);
     }else{
        tmpex = grpInd[nextH];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[nextH]); 
     }
   }

   // check the other neighbor
   prevO = search2(cid, atoms[cid].resNum, "O", 1);
   if(prevO>0){
     if(type==0){
        excl.push_back(prevO);
     }else{
        tmpex = grpInd[prevO];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[prevO]); 
     }
   }

   prevC = search2(cid, atoms[cid].resNum, "C", 1);
   if(prevC>0){
     if(type==0){
        excl.push_back(prevC);
     }else{
        tmpex = grpInd[prevC];
        if(find(excl.begin(), excl.end(), tmpex)==excl.end())
           excl.push_back(grpInd[prevC]);
     } 
   }

   prevN = search2(cid, atoms[cid].resNum, "N", 0);
   if(prevN>0){
      if(type==0){
         excl.push_back(prevN);
      }else{
         tmpex = grpInd[prevN];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[prevN]); 
      }
   }

   prevCa = search2(cid, atoms[cid].resNum, "CA", 0);
   if(prevCa>0){
      if(type==0){
         excl.push_back(prevCa);
      }else{
         tmpex = grpInd[prevCa];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[prevCa]);
      }
   }

   prevH = search2(cid, atoms[cid].resNum, "H", 0);
   if(prevH>0){
      if(type==0){
         excl.push_back(prevH);
      }else{
         tmpex = grpInd[prevH];
         if(find(excl.begin(), excl.end(), tmpex)==excl.end())
            excl.push_back(grpInd[prevH]);
      }
   }

   return excl;
}

void amideI::elst()
{
 
   int thisC, thisN;
   vector<int> this_exclude;

   rvec tmpEn, tmpEc;
   rvec Cxyz, Nxyz, tmpvec;

   float dc, dn, dj;

   fill_n(diag_w_el.begin(), nchrom, 0.0);

   // loop over all chromophores
   for(int i=0; i<nchrom; ++i){
      thisC = chrom_Clist[i];
      thisN = search2(thisC, atoms[thisC].resNum+1, "N", 1);
      
      this_exclude = get_excludes(thisC,0);
      // first pass; save excluded backbone atoms
      if(!save)
         exclude_list.push_back(this_exclude);

      setRvec(tmpEn,0.0);
      setRvec(tmpEc,0.0);      
      copyRvec(x[thisC],Cxyz);
      copyRvec(x[thisN],Nxyz);
      // loop over all atoms
      for(int j=0; j<natoms; ++j){
         // check if this atom belong to exclude list
         if(find(exclude_list[i].begin(), exclude_list[i].end(), j)!=exclude_list[i].end())
           continue; 

         // calc to dist to C atom
         addRvec(Cxyz, x[j], tmpvec, -1.0);
         pbc(tmpvec,box);
         dc = dist(tmpvec);
         
         // calc dist to N atom
         addRvec(Nxyz, x[j], tmpvec, -1.0);
         pbc(tmpvec,box);
         dn = dist(tmpvec);

         if(dc < constants::CN_dist_cutoff && dn < constants::CN_dist_cutoff){
            // atom j is within a cut-off distance to C and N
            if(atoms[j].charge){
               // calc Ec
               addRvec(Cxyz, x[j], tmpvec, -1);
               pbc(tmpvec,box);
               dj = dist(tmpvec);
               multRvec(tmpvec, atoms[j].charge/(dj*dj*dj));
               addRvec(tmpvec, tmpEc, 1);

               // calc En
               addRvec(Nxyz, x[j], tmpvec, -1);
               pbc(tmpvec,box);
               dj = dist(tmpvec);
               multRvec(tmpvec, atoms[j].charge/(dj*dj*dj));
               addRvec(tmpvec, tmpEn, 1);
            }
         }else if(dc < constants::CN_dist_cutoff){
            // just C atom is within a cut-off
            if(atoms[j].charge){
               // calc Ec
               addRvec(Cxyz, x[j], tmpvec, -1);
               pbc(tmpvec,box);
               dj = dist(tmpvec);
               multRvec(tmpvec, atoms[j].charge/(dj*dj*dj));
               addRvec(tmpvec, tmpEc, 1);
            }
         }else if(dn < constants::CN_dist_cutoff){
            // just N atom is within a cut-off
            if(atoms[j].charge){
               // calc En
               addRvec(Nxyz, x[j], tmpvec, -1);
               pbc(tmpvec,box);
               dj = dist(tmpvec);
               multRvec(tmpvec, atoms[j].charge/(dj*dj*dj));
               addRvec(tmpvec, tmpEn, 1);
            }
         }
      }
      // calc frequency
      diag_w_el[i] = map_el.getw01(dot(tmpEc, vecCO[i]), dot(tmpEn, vecCO[i]));
   }

   save = true;
}

void amideI::elstCG()
{
 
   int thisC, thisN;
   vector<int> this_exclude;

   rvec tmpEn, tmpEc;
   rvec Cxyz, Nxyz, tmpvec;
   rvec thisCcog, thisNcog;

   float dc, dn, dj;
   int kk;

   fill_n(diag_w_el.begin(), nchrom, 0.0);

   // loop over all chromophores
   for(int i=0; i<nchrom; ++i){
      thisC = chrom_Clist[i];
      thisN = search2(thisC, atoms[thisC].resNum+1, "N", 1);
      
      // here should be group excludes
      this_exclude = get_excludes(thisC,1);
      // first pass; save excluded backbone atoms
      if(!save)
         exclude_group_list.push_back(this_exclude);

      setRvec(tmpEn,0.0);
      setRvec(tmpEc,0.0);      

      copyRvec(cog[grpInd[thisC]],thisCcog);
      copyRvec(cog[grpInd[thisN]],thisNcog);

      copyRvec(x[thisC],Cxyz);
      copyRvec(x[thisN],Nxyz);
      // loop over all atoms
      for(int j=0; j<ncog; ++j){
         // check if this atom belong to exclude list
         if(find(exclude_group_list[i].begin(), exclude_group_list[i].end(), j)!=exclude_group_list[i].end())
            continue; 

         // calc to dist to C atom
         addRvec(thisCcog, cog[j], tmpvec, -1.0);
         pbc(tmpvec,box);
         dc = dist(tmpvec);
         
         // calc dist to N atom
         addRvec(thisNcog, cog[j], tmpvec, -1.0);
         pbc(tmpvec,box);
         dn = dist(tmpvec);

         if(dc < constants::CN_dist_cutoff && dn < constants::CN_dist_cutoff){
            // atom j is within a cut-off distance to C and N
            for(kk=chgSt[j]; kk<chgSt[j+1]; ++kk){
              if(atoms[kk].charge){
                 // calc Ec
                 addRvec(Cxyz, x[kk], tmpvec, -1);
                 pbc(tmpvec,box);
                 dj = dist(tmpvec);
                 multRvec(tmpvec, atoms[kk].charge/(dj*dj*dj));
                 addRvec(tmpvec, tmpEc, 1);

                 // calc En
                 addRvec(Nxyz, x[kk], tmpvec, -1);
                 pbc(tmpvec,box);
                 dj = dist(tmpvec);
                 multRvec(tmpvec, atoms[kk].charge/(dj*dj*dj));
                 addRvec(tmpvec, tmpEn, 1);
              }
            }
         }else if(dc < constants::CN_dist_cutoff){
            // just C atom is within a cut-off
            for(kk=chgSt[j]; kk<chgSt[j+1]; ++kk){
               if(atoms[kk].charge){
                  // calc Ec
                  addRvec(Cxyz, x[kk], tmpvec, -1);
                  pbc(tmpvec,box);
                  dj = dist(tmpvec);
                  multRvec(tmpvec, atoms[kk].charge/(dj*dj*dj));
                  addRvec(tmpvec, tmpEc, 1);
               }
            }
         }else if(dn < constants::CN_dist_cutoff){
            // just N atom is within a cut-off
            for(kk=chgSt[j]; kk<chgSt[j+1]; ++kk){
               if(atoms[kk].charge){
                  // calc En
                  addRvec(Nxyz, x[kk], tmpvec, -1);
                  pbc(tmpvec,box);
                  dj = dist(tmpvec);
                  multRvec(tmpvec, atoms[kk].charge/(dj*dj*dj));
                  addRvec(tmpvec, tmpEn, 1);
               }
            }
         }
      }
      // calc frequency
      diag_w_el[i] = map_el.getw01(dot(tmpEc, vecCO[i]), dot(tmpEn, vecCO[i]));
   }

   save = true;
}

void amideI::trDip()
{
//
// Calculate transition dipole moments
//
   int thisC, thisN;
   rvec Cxyz, Nxyz, n1, n2;
   fill_n(tdmuf.begin(), 3*nchrom, 0.0);

   for(int i=0; i<nchrom; ++i){
      setRvec(vecCO[i],0.0);
      setRvec(vecCN[i],0.0);
      setRvec(dip[i],0.0);

      thisC = chrom_Clist[i];
      // calc CO vector and its unit vector
      copyRvec(x[thisC],Cxyz);
      addRvec(x[thisC+1],Cxyz,vecCO[i],-1);
      pbc(vecCO[i],box);
      unitv(vecCO[i]);
      // calc CN vector and its unit vector
      thisN = search2(thisC, atoms[thisC].resNum+1, "N", 1);
      copyRvec(x[thisN],Nxyz);
      addRvec(Nxyz,Cxyz,vecCN[i],-1);
      pbc(vecCN[i],box);
      normalize(vecCN[i]);
      // calc TD
      cross(vecCO[i],vecCN[i],n1); //n1: normal to OCN plane
      cross(n1,vecCO[i],n2);    //n2: perpendicular to CO, in direction of N
      normalize(n2);

      copyRvec(n2,dip[i]);
      multRvec(dip[i],sin(constants::amideItdAngle));
      addRvec(vecCO[i],dip[i],-cos(constants::amideItdAngle)); // reverse CO vector

      normalize(dip[i]);
      multRvec(dip[i],constants::amideItdMag);
      //multRvec(dip[i],constants::amideItdMag);
      // copy into final array
      copyRRvec(&tdmuf[3*i],dip[i]);
      for(int jj=0; jj<3; ++jj)
         tdmuf[3*i+jj] *= constants::amideItdMag;

      // transition dipole location
      copyRvec(Cxyz,trdipL[i]);
      addRvec(vecCO[i],trdipL[i],constants::amideItd1);
      addRvec(vecCN[i],trdipL[i],constants::amideItd2); 
   }

}

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
      diag_w_nn[i]  = calc_N_NN_CS(thisC,thisRes,x,1);
      diag_w_nn[i] += calc_C_NN_CS(thisC,thisRes,x,1);
   }

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
  psi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);

  // calculate phi: C-N-Ca-C
  i1=search2(atomI,resI-1,"C",1);
  i2=search2(atomI,resI,"N",0);
  i3=search2(atomI,resI,"CA",0);
  i4=atomI;
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

float amideI::calc_N_NN_CS(const int thisC, const int thisRes, const rvec *x, int mode)
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

   // calculate phi: C-N-Ca-C
   i1=search2(thisC,thisRes,"C",1);
   i2=search2(thisC,thisRes,"N",0);
   i3=search2(thisC,thisRes,"CA",0);
   i4=thisC; //atomI;
   phi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);

   if(mode==1){
      return map_nn.getNNshift(phi,psi,"Nterm");
   }else{
      return map_nnc.getCoupling(phi,psi);
   }

   //return map_nn.getNNshift(phi,psi,"Nterm");
}


float amideI::calc_C_NN_CS(const int thisC, const int thisRes, const rvec *x, int mode)
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
   //cout << " C shift atoms psi " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

   // calculate phi: C-N-Ca-C
   i1=thisC; //search2(thisC,resI-1,"C",1);
   i2=search2(thisC,thisRes+1,"N",1);
   i3=search2(thisC,thisRes+1,"CA",1);
   i4=search2(thisC,thisRes+1,"C",1); //atomI;
   phi = calcDihedral(x[i1],x[i2],x[i3],x[i4]);

   if(mode==1){
      return map_nn.getNNshift(phi,psi,"Cterm");
   }else{
      return map_nnc.getCoupling(phi,psi);
   }
}

void amideI::updateExd()
{
//
// set up the diagonal elements of excitonic Hamiltonian here
//
   fill_n(hf.begin(), nchrom2, 0.0);
   for(int ii=0; ii<nchrom; ++ii)
      hf[ii*nchrom+ii] = diag_w_nn[ii] + diag_w_el[ii] + isoShift;
   
}

amideI::~amideI() {
   delete [] vecCO;
   delete [] vecCN;
   delete [] trdipL;
   delete [] dip;
   delete [] cog; 
}


void amideI::writeH()
{
  int jj = 0;
  int pr_size;
  pr_size = nchrom;
  for(int ii=0; ii<nchrom; ++ii){
     jj = ii*nchrom+ii;
     houtfile.write(reinterpret_cast<char*>(&hf[jj]), pr_size*sizeof(float));
     pr_size -= 1;
  }

  //for(int i=0; i<nchrom; ++i)
  //   for(int j=i; j<nchrom; ++j)
  //      cout << hf[i*nchrom+j] << " ";
  //cout << endl;
}

void amideI::writeJ()
{ jobfile.write(reinterpret_cast<char*>(&nchrom), sizeof(int)); }

void amideI::writeD()
{ doutfile.write(reinterpret_cast<char*>(&tdmuf[0]), (3*nchrom)*sizeof(float)); 

  //for(int i=0; i<(3*nchrom); ++i)
  //    cout << tdmuf[i] << " ";
  //cout << endl;
}


