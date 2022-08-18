#include <cstdio>
#include <string>
#include <iostream>
#include <chrono>
#include "input.h"
#include "amideI.h"
#include "../util/util.h"

using namespace std;

int main(int argc, char ** argv){

   printv();
   printf("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Amide I module <<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
   tstamp("\n** Program starts @ ");
   uhstName();

   // read input parameters
   Input input(argc, argv);   

   std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

   // create and run simulation
   amideI S (input.getGroFile(), input.getTrajFile(), input.getItpfs(),
             input.getTopFile(), input.getJobType());

   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

   timet(start, end);
   tstamp("\n** Program finished @ ");   
   printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Done. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");

   return EXIT_SUCCESS;
}
