#include "util.h"

using namespace std;

void printv()
{
   printf("\n\n");
   printf("*************************************************************************\n");
   printf("*                                                                       *\n");
   printf("*        MultiSpec: a Package for Simulating Vibrational Spectra        *\n");
   printf("*                       of Condensed-phase Systems                      *\n");
   printf("*                                                                       *\n");
   printf("*                              Version 3.1                              *\n");
   printf("*                              March  2024                              *\n");
   printf("*              http://github.com/kananenka-group/MultiSpec              *\n");
   printf("*                                                                       *\n");
   printf("*             Copyright (c) 2021 - 2024 Alexei A. Kananenka             *\n");
   printf("*                                                                       *\n");
   printf("*************************************************************************");

}

void fileReadErr(string _fn_)
{
   printf("Error! Reading file %s file. Probably reached EOF.\n",_fn_.c_str());
   exit(EXIT_FAILURE);
}

void remove_leading_trailing_whitespace(string &str)
{ 
  str.erase(str.begin(), find_if(str.begin(), str.end(), bind1st(std::not_equal_to<char>(), ' '))); 
  str.erase(find_if(str.rbegin(), str.rend(), bind1st(not_equal_to<char>(), ' ')).base(), str.end());
}

string str_toupper(string s) {
    transform(s.begin(), s.end(), s.begin(), 
              [](unsigned char c){ return toupper(c); } 
              );
    return s;
}

float switchf(float ez)
{
//
// Calculate switching function here.
// See J. Chem. Phys. 135, 044701 (2011)
// 
   float rc, rc2, rc3;
   float fz1, fz2;

   rc  = constants::SWITCHF_CUT*constants::A0;
   rc2 = rc*rc;
   rc3 = rc2*rc;
   fz1 = 0;
   fz2 = 0;

   fz1 = 1;
   if(ez <= rc && ez >= -rc)
     fz1 = (2.0*rc3 + 3.0*rc2*ez - ez*ez*ez)/(4.0*rc3);
   
   if(ez <= -rc) fz1 = 0;
   
   fz2 = -1;
   ez = -ez;
   if(ez <= rc && ez >= -rc)
     fz2 = -(2.0*rc3 + 3.0*rc2*ez - ez*ez*ez)/(4.0*rc3);
   
   if(ez <= -rc)  fz2 = 0;
   
   return (fz1 + fz2);
}

void removeFile(vector<string> filenames)
{
   int result;
   printf("\n** Old temporary files found and deleted : "); 
   for(unsigned int ii=0; ii<filenames.size(); ++ii){
      result = remove(filenames[ii].c_str());
      if(result == 0 ) //{
         printf(" %s ",filenames[ii].c_str()); 
   }
   printf("\n");
}

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
