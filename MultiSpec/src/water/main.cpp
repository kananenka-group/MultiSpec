#include <cstdio>
#include <string>
#include <iostream>
#include <chrono>
#include "input.h"
#include "water.h"
#include "../time/timeSt.h"

using namespace std;

int main(int argc, char ** argv){

   printf("\n\n>>>>>>>>>>>>>>>>>>>> Water module <<<<<<<<<<<<<<<<<<<< \n");

   // read input parameters
   Input input(argc, argv);   

   std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

   // create and run simulation
   water S (input.getWaterModel(), input.getChromType(), input.getNFrames(), 
            input.getOHMap(), input.getJobType(), input.getTrajFile(), 
            input.getGroFile(), input.getChgFile(), input.ifIR(), input.ifRaman(),
            input.ifSFG());

   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   
   timet(start, end);
   printf("\n>>>>>>>>>>>>>>>>>>>> Done. <<<<<<<<<<<<<<<<<<<<\n\n");

   return EXIT_SUCCESS;
}
