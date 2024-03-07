#include "traj.h"

using namespace std;

Traj::Traj(const char *xtcfile)
{
  printf("\n** Checking trajectory file: %s. **\n",xtcfile);

  //get number of atoms
  read_xtc_natoms((char*)xtcfile,&natoms);

  //open xtc file
  trj = xdrfile_open(xtcfile,"r");
  if (trj == NULL) {
    printf("ERROR: trajectory file %s cannot be read.\n",xtcfile);
    exit(EXIT_FAILURE);
  } 

  //initialize vars
  nT = 0;
  x  = new rvec[natoms];
}

Traj::~Traj()
{
  delete[] x;
  xdrfile_close(trj);
}

//read next step of trajectory
int Traj::next(const bool convertFlag)
{
  //read one timestep, returns 0 if success
  int stat=read_xtc(trj,natoms,&step,&t,boxT,x,&prec);

  //return if finished reading file
  if (stat) return stat;

  if (convertFlag) convert();

  nT++;

  return stat;
}

void Traj::convert()
{
  for (int jj=0; jj<natoms; jj++)
    for (int kk=0; kk<DIM; kk++)
      x[jj][kk]*=A0INV;

  for (int ii=0; ii<DIM; ii++)
    box[ii]=boxT[ii][ii]*A0INV;
}


void Traj::getBox(rvec &out) const {
  for (int kk=0; kk<DIM; kk++)
    out[kk]=box[kk];
}

