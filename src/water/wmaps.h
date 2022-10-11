#ifndef WMAPS_H
#define WMAPS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <functional>
#include "mapsTs.h"
#include "../vec/vecManip.h"
#include "../util/util.h"

using namespace std;

class WmapS{

public:
  WmapS(string, string, bool);
  ~WmapS() {};

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
  float getcnn(const float, const float, const float,const float, const float, const float);

private:

  eFmapS emap; 
  string map_name;
  string map_name_inp;
  string ctype;

  bool intrac;

};


#endif
