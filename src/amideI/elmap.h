#ifndef ELMAP_H
#define ELMAP_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "../util/util.h"

using namespace std;

class ELmap{

public:
  ELmap(string, string);
  ~ELmap() {};

  float getw01(const float ec, const float en);

private:

  string map_name;
  string map_name_inp;
  string res_map_file;

  int mapi;

  void checkResMaps();
  
};

#endif
