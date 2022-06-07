#!/bin/bash -l

stat=10
frames=$((stat*350))

srcdir='/work/akanane/sw/MultiSpec/src'
tstdir='/work/akanane/sw/MultiSpec/tests'

$srcdir/water/water_gen --xtc $tstdir/traj.xtc  --gro_file $tstdir/confout.gro --charge_file $tstdir/q.dat --OHmap gruenbaum_tip4p_2013_OH  --IR 0 --Raman 1  --nframes $frames --chrom_type ws  --spec_type pure  --water_model tip4p

$srcdir/exciton/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --P Polarizability.bin  --IR 0 --Raman 1  --nframes $frames  --nchrom 1000  --T1 0.26  --navg $stat  --tsep 2.0

