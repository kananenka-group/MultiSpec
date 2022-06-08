#ifndef INPUT_H
#define INPUT_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class Input{

public:
   Input(int argc, char ** argv);
  ~Input() {};

  string getTrajFile()   const { return xtc_file; }
  string getOHMap()      const { return OHmap_name; }
  string getOHMapA()     const { return ODmap_name; }
  string getChromType()  const { return chrom_type; }
  string getWaterModel() const { return water_model; }
  string getGroFile()    const { return gro_file; }
  string getJobType()    const { return job_type; }
  string getChgFile()    const { return chg_file; }
  int getNFrames()       const { return nframes;  }
  int getND2O()          const { return nd2o; }
  bool ifIR()            const { return ir; }
  bool ifRaman()         const { return raman; }
  bool ifSFG()           const { return sfg; }

private:
  int nframes;
  int nd2o;

  bool ir    = true; 
  bool raman = false; 
  bool sfg   = false;

  string xtc_file    = "traj.xtc";
  string chrom_type  = "OH_stretch";
  string OHmap_name  = "gruenbaum_tip4p_2013_OH";
  string ODmap_name  = "gruenbaum_tip4p_2013_OD";
  string job_type;
  string water_model = "tip4p";;
  string gro_file    = "confout.gro";
  string chg_file    = "charges.dat";

};

#endif
