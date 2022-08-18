#ifndef INPUT_H
#define INPUT_H

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
  vector<string> getItpfs() const { return itpfs; }

private:
  int nframes = 1;
  int startframe = 1;

  bool ir    = false; 
  bool excHam = false;

  vector<string> itpfs;

  string xtc_file    = "traj.xtc";
  string job_type    = "full";
  string gro_file    = "confout.gro";
  string top_file    = "topol.top";

  ofstream inpfile;

};

#endif
