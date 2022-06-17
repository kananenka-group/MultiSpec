#ifndef WMAPB_H
#define WMAPB_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <functional>
#include "mapsTb.h"
#include "../vec/vecManip.h"

using namespace std;

class WmapB{

public:
  WmapB(string, string);
  ~WmapB() {};

  float getw01E(const float);
  float getw12E(const float);
  float getw02E(const float);
  float getw01E_HOH(const float);
  float getw12E_HOH(const float);
  float getw01E_DOD(const float);
  float getw12E_DOD(const float);
  float getw01E_HOD(const float);
  float getw12E_HOD(const float);

private:

  eFmapB emap; 
  string ctype;

};


#endif
