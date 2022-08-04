#!/bin/bash -l

stat=1
frames=$((stat*350))

srcdir='/work/akanane/sw/MultiSpec/src/build'
tstdir='/work/akanane/users/akanane/Raman_Water/traj/E3B2'

$srcdir/water_gen --xtc $tstdir/traj.xtc  --gro_file $tstdir/confout.gro --atoms_file q.inp --stretch_map gruenbaum_2013_tip4p  --IR 1 --Raman 0 --nframes $frames   --spec_type wsOH --water_model e3b2

$srcdir/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --D Dipole.bin  --IR 1 --Raman 0 --nframes $frames  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2 --inh 1

