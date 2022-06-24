#include "util.h"

using namespace std;

void tstamp(std::string message)
{
   string timem = currentDateTime();
   cout << message << " " << timem << endl;
}

void timet(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
   double timet = (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
   timet /= 1000;
   int hrs = floor(timet/3600.0);
   int min = floor(timet/60.0) - hrs*60.0;
   double sec = timet - hrs*3600 - min*60;
   printf("\n** Total time: %d hours  %d minutes %d seconds.\n",hrs,min,(int)round(sec));
}

const std::string currentDateTime()
{
   time_t     now = time(0);
   struct tm  tstruct;
   char       buf[80];
   tstruct = *localtime(&now);
   strftime(buf, sizeof(buf), "%Y/%m/%d %X", &tstruct);
   return buf;
}

void uhstName()
{
   char hostname[HOST_NAME_MAX];
   char username[LOGIN_NAME_MAX];
   gethostname(hostname, HOST_NAME_MAX);
   getlogin_r(username, LOGIN_NAME_MAX);
   //printf("   Username: %s \n",username);
   printf("   Hostname: %s \n",hostname);
}
