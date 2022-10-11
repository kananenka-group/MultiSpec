#include "input.h"

namespace po = boost::program_options;
using namespace std;

Input::Input(int argc, char ** argv){

   // parse input options
   try {

       po::options_description desc("Allowed options");
       desc.add_options()
       ("help", "help message")
       ("xtc",            po::value<string>(&xtc_file),         "gromacs *.xtc file")
       ("IR",             po::value<bool>(&ir),                 "Calculate linear IR spectrum")
       ("nframes",        po::value<int>(&nframes),             "Number of frames to read")
       ("spec_type",      po::value<string>(&job_type),         "Type of spectrum")
       ("isotope_labels", po::value<vector<string> >(&isolabels)->multitoken(), "isotope labeled residues") 
       ("gro_file",       po::value<string>(&gro_file),         "gromacs *.gro file")
       ("start",          po::value<int>(&startframe),          "starting frame to read from xtc file")
       ("exc_ham",        po::value<bool>(&excHam),             "print diagonal frequencies and couplings")
       // https://stackoverflow.com/questions/42428671/boost-program-options-multiple-values-for-an-option
       // multiple values for boost option
       ("itp_files",      po::value<vector<string> >(&itpfs)->multitoken(),   "itp files (provide all of them)")
       ("top_file",       po::value<string>(&top_file),         "gromacs *.top file")
       ("nn_map",         po::value<string>(&nn_map),           "nearest-neighbor frequency map")
       ("nnc_map",        po::value<string>(&nnc_map),           "nearest-neighbor coupling map")
       ("el_map",         po::value<string>(&el_map),           "electrostatic backbone map")
       ("isotope_shift",  po::value<float>(&isoShift),          "isotope frequency shift")
       ("residues_maps",  po::value<string>(&res_map_file),     "file containing maps for specific residues")
       ("exc_ham",        po::value<bool>(&excHam),       "print diagonal frequencies, inter- and intramolecular couplings")
       //("cut_off_elst",   po::value<string>(&cutoffelst),       "cut-off for calculating electrostatics")
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

    // print input parameters
    inpfile.open("amideI.inp");
    inpfile << "# Input parameters of amideI module " << endl;
    inpfile << "           xtc = " << xtc_file << endl;
    inpfile << "            IR = " << ir << endl;
    inpfile << "       nframes = " << nframes << endl;
    inpfile << "     spec_type = " << job_type << endl;
    inpfile << "      gro_file = " << gro_file << endl;
    inpfile << "         start = " << startframe << endl;
    inpfile << "       exc_ham = " << excHam << endl;
    inpfile << "      top_file = " << top_file << endl;
    inpfile << "        nn_map = " << nn_map << endl;
    inpfile << "       nnc_map = " << nnc_map << endl;
    inpfile << "        el_map = " << el_map << endl;
    inpfile << "     itp_files = ";
    inpfile << " isotope_shift = " << isoShift << endl;
    //inpfile << "  cut_off_elst = " << cutoffelst << endl;
    for(uint i=0; i<itpfs.size(); ++i)
       inpfile << itpfs[i] << " ";
    inpfile << endl;
    if(job_type=="iso"){
       inpfile << "isotope_labels = ";
       for(uint i=0; i<isolabels.size(); ++i)
          inpfile << isolabels[i] << " ";
       inpfile << endl;
    }
    if(!res_map_file.empty())
       inpfile << " residues_maps = " << res_map_file << endl;
    inpfile.close();

}

