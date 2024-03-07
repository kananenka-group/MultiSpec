#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <boost/program_options.hpp>

using namespace std;

class Input{

public:
   Input(int argc, char ** argv);
  ~Input() {};

  string getTrajFile()   const { return xtc_file; }
  string getwSMap()      const { return stretch_map_name; }
  string getwBMap()      const { return bend_map_name; }
  string getWaterModel() const { return water_model; }
  string getGroFile()    const { return gro_file; }
  string getJobType()    const { return job_type; }
  string getAmsFile()    const { return ams_file; }
  int getNFrames()       const { return nframes;  }
  int getND2O()          const { return nd2o; }
  bool ifIR()            const { return ir; }
  bool ifRaman()         const { return raman; }
  bool ifSFG()           const { return sfg; }
  bool ifDODv()          const { return dodov; }
  bool getIntraC()       const { return intrac; }
  bool getInermcs()      const { return intercs;}
  float getFc()          const { return fc; }
  float getTDSFG()       const { return trdipSFG; }
  int getStartFrame()    const { return startframe; }
  bool getExcHam()       const { return excHam; }
  float getImcut()       const { return imcut; }       
  float getwmin()        const { return wmin; }
  float getwmax()        const { return wmax; }
  int getNbins()       const { return nbins; }

private:
  int nframes = 1;
  int nd2o = 0;
  int startframe = 1;

  bool ir    = false; 
  bool raman = false; 
  bool sfg   = false;
 
  float fc = 25.0;  // OH stretch-HOH bend overtone Fermi coupling
  float trdipSFG=-10;
  float imcut=1.0e6;
  float wmin = 1000;
  float wmax = 5000;
  int nbins = 200;

  bool dodov = false;
  bool intrac= true;
  bool intercs=true;

  bool excHam = false;

  string xtc_file    = "traj.xtc";
  string stretch_map_name; 
  string bend_map_name;
  string job_type    = "wsOH";
  string water_model = "tip4p";
  string gro_file    = "confout.gro";
  string ams_file    = "atoms.txt";

  ofstream inpfile;

};

#endif
