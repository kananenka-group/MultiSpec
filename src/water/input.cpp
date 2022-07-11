#include "input.h"

namespace po = boost::program_options;
using namespace std;

Input::Input(int argc, char ** argv){

   // parse input options
   try {

       po::options_description desc("Allowed options");
       desc.add_options()
       ("help", "help message")
       ("xtc",         po::value<string>(&xtc_file),   "xtc file")
       ("stretch_map", po::value<string>(&stretch_map_name), "stretch map")
       ("bend_map",    po::value<string>(&bend_map_name), "bend map")
       ("water_model", po::value<string>(&water_model),"water model")
       ("IR",          po::value<bool>(&ir),           "Calculate linear IR spectrum")
       ("Raman",       po::value<bool>(&raman),        "Calculate Raman spectrum") 
       ("SFG",         po::value<bool>(&sfg),          "Calculate SFG spectrum")
       ("nframes",     po::value<int>(&nframes),       "Number of frames to read")
       ("spec_type",   po::value<string>(&job_type),   "Type of spectrum")
       ("gro_file",    po::value<string>(&gro_file),   "gromacs file")
       ("atoms_file",  po::value<string>(&ams_file),   "file with atomic information: charges, masses")
       ("D2O",         po::value<int>(&nd2o),          "number of D2O molecules")
       ("DOD_overtone",po::value<bool>(&dodov),        "turn on/off DOD overtone in iso mix simulations")
       ("Fc",          po::value<float>(&fc),          "stretch-bend Fermi coupling")
       ("trdipSFG",    po::value<float>(&trdipSFG),    "distance from O atom and transition dipole for SFG")
       ("intrac",      po::value<bool>(&intrac),       "intramolecular OH stretch coupling")
       ("intercOH",    po::value<bool>(&intercs),      "OH stretch intermolecular coupling")
       ;

       po::variables_map vm;
       po::store(po::parse_command_line(argc, argv, desc), vm);
       po::notify(vm);

      if (vm.count("help")) {
         cout << " Help message\n" << endl;
         exit(EXIT_FAILURE);
       }


    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        exit(EXIT_FAILURE);
    }

}

