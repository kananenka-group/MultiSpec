#include <cstdio>
#include <string>
#include <chrono>
#include "input.h"
#include "exc.h"
#include "../time/timeSt.h"

using namespace std;

int main(int argc, char ** argv){

   printf("\n\n>>>>>>>>>>>>>>>>>>>> Exciton module <<<<<<<<<<<<<<<<<<<< \n");

   // read input parameters
   Input input(argc, argv);   
 
   // run simulation
   Exc S (input.getHfile(), input.getDfile(), input.getNchrom(),
          input.getNfrm(), input.getdt(), input.gettc(), 
          input.getTrlx(), input.getNavg(), input.getTsep(),
          input.ifIR(), input.ifRaman(), input.getWavg());

   std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

   S.run();

   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

   // print running time
   timet(start, end);

   printf("\n>>>>>>>>>>>>>>>>>>>> Done. <<<<<<<<<<<<<<<<<<<<\n\n");

   return EXIT_SUCCESS;
}
