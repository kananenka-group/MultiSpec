#!/bin/bash -l

stat=1
frames=$((stat*350))

srcdir='/work/akanane/sw/MultiSpec/src'
tstdir='/work/akanane/users/akanane/SFG/sfg_test'

$srcdir/water/water_gen --xtc $tstdir/sfg_test.xtc  --gro_file $tstdir/confout.gro --atoms_file q.inp --stretch_map gruenbaum_2013_tip4p  --SFG 1 --nframes $frames --spec_type wsOH --water_model tip4p 

$srcdir/exciton/excp --dt 0.01  --tc 1.2  --H Hamiltonian.bin  --D Dipole.bin  --P Polarizability.bin --SFG 1 --nframes $frames  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2

