#ifndef TRAJ_H
#define TRAJ_H

#include <xdrfile_xtc.h>
#include <xdrfile.h>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <cmath>

#define A0INV 18.89726125    //inverse bohr radius (nm/a0)
#define PI 3.14159265359

using namespace std;
class Traj
{
public:
  Traj(const char *xtcfile);
  ~Traj();
  int next(const bool convertFlag=true);

  int   getNatoms() const {return natoms;};
  int   getNT()     const {return nT;};
  float getT()      const {return t;};

  const rvec* getCoords() const { return x; };
  void getBox(rvec &box) const;

private:
  XDRFILE *trj;   //pointer to file
  int natoms;     //number of atoms
  int nT;         //number of timesteps processed

  int step;       //integer timestep
  float t;        //time in ps (units correct?)
  float prec;     //precision
  matrix boxT;    //simulation box in raw output
  rvec box;       //simulation box
  rvec *x;        //atom positions

  float dt;       //timestep between current step and last

  void convert();
};

#endif
