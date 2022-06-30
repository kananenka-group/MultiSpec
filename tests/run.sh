#!/bin/bash -l

stat=1
frames=$((stat*350))

srcdir='/work/akanane/sw/MultiSpec/src'
tstdir='/work/akanane/users/akanane/Raman_Water/traj'

$srcdir/water/water_gen --xtc $tstdir/traj.xtc  --gro_file $tstdir/confout.gro --atoms_file q.inp --stretch_map gruenbaum_2013_tip4p  --IR 1 --Raman 1 --nframes $frames   --spec_type wswbD2O --water_model tip4p 

$srcdir/exciton/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --D Dipole.bin  --IR 1 --Raman 1 --nframes $frames  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2

