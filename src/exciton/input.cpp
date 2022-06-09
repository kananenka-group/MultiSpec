#include "input.h"
#include <exception>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

Input::Input(int argc, char ** argv){

   // parse input options
   try {

       po::options_description desc("Allowed options");
       desc.add_options()
       ("help", "help message")
       ("dt",     po::value<double>(&dt),              "Time step [ps]")
       ("tc",     po::value<double>(&tc),              "Propagation time [ps]")
       ("tsep",   po::value<double>(&tsep),            "Separation time [ps]")
       ("T1",     po::value<double>(&trlx),            "Relaxation time [ps]")
       ("w_avg",  po::value<double>(&w_avg),           "Average frequency [cm-1]")
       ("H",      po::value<std::string>(&Hfile),      "Hamiltonian file")
       ("D",      po::value<std::string>(&Dfile),      "Dipole file")
       ("P",      po::value<std::string>(&Pfile),      "Polarizability file")
       ("nchrom", po::value<int>(&nchrom),             "Number of chromophores")
       ("nframes",po::value<int>(&nframe),             "Number of frames to read")
       ("navg",   po::value<int>(&navg),               "Trajectories for statistical averaging")
       ("IR",     po::value<bool>(&ir),                "Calculate linear IR spectrum")
       ("Raman",  po::value<bool>(&raman),             "Calculate Raman spectrum") 
       ("SFG",    po::value<bool>(&sfg),               "Calculate SFG spectrum")
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


