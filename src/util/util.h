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
#include <algorithm>
#include <limits.h>
#include "../const/const.h"

using namespace std;

void remove_leading_trailing_whitespace(string &str);
void printv();
void tstamp(std::string message);
void timet(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end);
const std::string currentDateTime();
void uhstName();
void removeFile(vector<string> filenames);
float switchf(float ez);
string str_toupper(string s);

#endif
