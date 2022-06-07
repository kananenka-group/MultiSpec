#include "nummath.h"

using namespace std;

double simpsonInt(double a, double b,
                    int n, func_type f)
{
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
       sum_odds += f(std::fma(i,h,a));
    }
    double sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
       sum_evens += f(std::fma(i,h,a);
    }

    return (std::fma(2,sum_evens,f(a)) + std::fma(4,sum_odds,f(b))) * h / 3;
}
