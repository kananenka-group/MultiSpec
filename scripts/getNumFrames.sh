#!/bin/bash

######################################################
#                                                    #
# Gromacs script for getting the total number        #
# of frames in a trajectory file                     #
#                                                    #
######################################################
FILE="../tests/trajL.xtc"
vpkg_require gromacs/2018.1
gmx check -f "$FILE"
