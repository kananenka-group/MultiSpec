#ifndef TIME_H
#define TIME_H

#include <string>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include <chrono>
#include <unistd.h>
#include <limits.h>

void tstamp(std::string message);
void timet(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end);
const std::string currentDateTime();
void uhstName();

#endif
