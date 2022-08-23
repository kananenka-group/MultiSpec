#ifndef NNMAP_H
#define NNMAP_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

class NNmap{

public:
  NNmap(string);
  ~NNmap() {};

  //vector<float> getNtermShifts() { return NtermShift; }
  //vector<float> getCtermShifts() { return CtermShift; }
  float getdT()  { return dTheta; }
  int getnT()    { return nTheta; }
  float getNNshift(const float &x, const float &y, const string termshift);

private:

  string map_name;
  
  vector<float> NtermShift;
  vector<float> CtermShift;

  float dTheta;
  int nTheta;


};

#endif
