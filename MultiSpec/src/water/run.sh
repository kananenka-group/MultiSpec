#!/bin/bash

./water_gen --xtc /work/akanane/sw/spec/tests/traj.xtc \
            --gro_file /work/akanane/sw/spec/tests/confout.gro \
            --charge_file /work/akanane/sw/spec/tests/q.dat \
            --map gruenbaum_tip4p_2013_OH \
            --IR 1 \
            --nframes 1000 \
            --chrom_type ws \
            --spec_type pure \
            --water_model tip4p 
