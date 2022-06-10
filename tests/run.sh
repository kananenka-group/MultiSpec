#!/bin/bash -l

stat=1
frames=$((stat*350))

srcdir='/work/akanane/sw/MultiSpec/src'
tstdir='/work/akanane/sw/MultiSpec/tests'

$srcdir/water/water_gen --xtc $tstdir/traj.xtc  --gro_file $tstdir/confout.gro --charge_file $tstdir/q.inp --stretch_map gruenbaum_tip4p_2013  --IR 1 --Raman 0 --nframes $frames --chrom_type ws  --spec_type wsOH  --water_model tip4p 

$srcdir/exciton/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --D Dipole.bin  --IR 1 --Raman 0 --nframes $frames  --nchrom 1000  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2

