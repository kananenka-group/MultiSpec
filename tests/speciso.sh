#!/bin/bash -l

stat=100
frames=$((stat*350))

# fermi: true/false
fermi=true
Raman=1
IR=1
fermicoupling=25
D2O=100

srcdir='/work/akanane/sw/MultiSpec/src'
tstdir='/work/akanane/sw/MultiSpec/tests'
traj='trajL.xtc'

if $fermi ; then
   specType='wswbiso'
   $srcdir/water/water_gen --xtc $tstdir/$traj  --gro_file $tstdir/confout.gro --charge_file $tstdir/q.inp --stretch_map gruenbaum_2013_tip4p  --IR $IR --Raman $Raman --nframes $frames   --spec_type $specType --Fc $fermicoupling  --water_model tip4p --bend_map ni_2015_kananenka_2019_tip4p --D2O=$D2O --DOD_overtone 0
   $srcdir/exciton/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --D Dipole.bin  --P Polarizability.bin --IR $IR --Raman $Raman --nframes $frames  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2
else
   specType='wsiso'
   $srcdir/water/water_gen --xtc $tstdir/$traj  --gro_file $tstdir/confout.gro --charge_file $tstdir/q.inp --stretch_map gruenbaum_2013_tip4p  --IR $IR --Raman $Raman --nframes $frames  --spec_type $specType --water_model tip4p --D2O=$D2O
   $srcdir/exciton/excp --dt 0.01  --tc 1.5  --H Hamiltonian.bin  --D Dipole.bin  --IR $IR --P Polarizability.bin --Raman $Raman --nframes $frames  --T1 0.26  --navg $stat  --tsep 2.0 --w_avg 3415.2
fi

