#include <cstdio>
#include <string>
#include <iostream>
#include <chrono>
#include "input.h"
#include "water.h"
#include "../util/util.h"

using namespace std;

int main(int argc, char ** argv){

   printv();
   printf("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Water module <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
   tstamp("\n** Program starts @ ");
   uhstName();

   // read input parameters
   Input input(argc, argv);   

   std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

   // create and run simulation
   water S (input.getWaterModel(), input.getwSMap(), input.getwBMap(),
            input.getJobType(), input.getTrajFile(), input.getGroFile(),
            input.getAmsFile(), input.getNFrames(), input.getND2O(),
            input.getStartFrame(), input.ifIR(), input.ifRaman(), 
            input.ifSFG(), input.ifDODv(), input.getExcHam(),
            input.getIntraC(), input.getInermcs(), input.getTDSFG(), 
            input.getFc());

   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

   timet(start, end);
   tstamp("\n** Program finished @ ");   
   printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Done. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");

   return EXIT_SUCCESS;
}
