#ifndef TIME_H
#define TIME_H

#include <chrono>

void timet(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
   double timet = (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
   timet /= 1000;
   int hrs = floor(timet/3600.0);
   int min = floor(timet/60.0) - hrs*60.0;
   double sec = timet - hrs*3600 - min*60;
   printf("\n Total time: %d hours  %d minutes %3.1f seconds.\n",hrs,min,sec);
}

#endif
