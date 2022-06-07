#ifndef WMAPS_H
#define WMAPS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "mapsT.h"
#include "../vec/vecManip.h"

using namespace std;

class Wmap{

public:
  Wmap(string);
  ~Wmap() {};

  float getw01E(const float); 
  float getm01E(const float); 
  float getx01E(const float);
  float getp01E(const float);
  float getcnn(const float, const float); 

private:

  eFmap emap; 

};


#endif
