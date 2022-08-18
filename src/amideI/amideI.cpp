#include "amideI.h"

amideI::amideI(string gro_file, string traj_file, vector<string> itp_files,
               string top_file, string spec_type) :
               s(gro_file, itp_files, top_file), traj_file(traj_file), spec_type(spec_type)
{




}

amideI::~amideI() { }
