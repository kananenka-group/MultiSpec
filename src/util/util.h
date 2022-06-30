#ifndef TIME_H
#define TIME_H

#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>
#include <unistd.h>
#include <limits.h>
#include "../const/const.h"

using namespace std;

void tstamp(std::string message);
void timet(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end);
const std::string currentDateTime();
void uhstName();
void removeFile(vector<string> filenames);
float switchf(float ez);

#endif
