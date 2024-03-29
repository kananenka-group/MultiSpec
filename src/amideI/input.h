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

  string getTrajFile()      const { return xtc_file; }
  string getGroFile()       const { return gro_file; }
  string getTopFile()       const { return top_file; }
  string getJobType()       const { return job_type; }
  int getNFrames()          const { return nframes;  }
  bool ifIR()               const { return ir; }
  int getStartFrame()       const { return startframe; }
  bool getExcHam()          const { return excHam; }
  bool getNISE()             const { return nise; }
  vector<string> getItpfs() const { return itpfs; }
  vector<string> getIsoL()  const { return isolabels; }
  string getNNmap()         const { return nn_map; }
  string getNNCmap()        const { return nnc_map; }
  string getELmap()         const { return el_map; }
  string getResMapFile()    const { return res_map_file;}
  float getIsoShift()       const { return isoShift; }

private:
  int nframes = 1;
  int startframe = 1;

  bool ir    = false; 
  bool excHam = false;
  bool nise = false;

  // 13C=18O isotope shift in cm-1
  // taken from JACS, 134, 19118-19128 (2012)
  float isoShift = -66.0;

  vector<string> itpfs;
  vector<string> isolabels;

  string xtc_file    = "traj.xtc";
  string job_type    = "full";
  string gro_file    = "confout.gro";
  string top_file    = "topol.top";
  string nn_map      = "jansen_2006";
  string nnc_map     = "jansen_2006";
  string el_map      = "wang_2011";
  string res_map_file= "";
  // implement different cut-offs for electrostatics later
  //string cutoffelst  = "cg";

  ofstream inpfile;

};

#endif
