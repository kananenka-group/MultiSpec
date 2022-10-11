#ifndef NNCMAP_H
#define NNCMAP_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "../util/util.h"

using namespace std;

class NNCmap{

public:
  NNCmap(string);
  ~NNCmap() {};

  float getCoupling(const float &,const float &);

private:

  string map_name;
  string map_name_inp;

  int nTheta;
  float dTheta;

  vector<float> coupling;

  
};

#endif
