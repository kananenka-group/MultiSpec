#ifndef WMAPS_H
#define WMAPS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <functional>
#include "mapsT.h"
#include "../vec/vecManip.h"

using namespace std;

class Wmap{

public:
  Wmap(string, string);
  ~Wmap() {};

  float getw01E(const float);
  float getm01E(const float);
  float getx01E(const float);
  float getp01E(const float);

  float getw01E_OH(const float); 
  float getm01E_OH(const float); 
  float getx01E_OH(const float);
  float getp01E_OH(const float);

  float getw01E_OD(const float); 
  float getm01E_OD(const float); 
  float getx01E_OD(const float);
  float getp01E_OD(const float);
  float getcnn(const float, const float); 

private:

  eFmap emap; 

  string ctype;

};


#endif
