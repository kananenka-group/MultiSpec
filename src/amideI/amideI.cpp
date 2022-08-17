#include "amideI.h"

amideI::amideI(string gro_file, string traj_file, vector<string> itp_files) :
               s(gro_file, itp_files), traj_file(traj_file)
{




}

amideI::~amideI() { }
