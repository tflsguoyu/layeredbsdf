#pragma once
#include <cmath>

class Statistics
{
public:
    Statistics(const Statistics &stats);
    Statistics()                    { reset(); }

    void reset()                    { n = 0; mean = M2 = 0.0; }
    void push(double x);
    void push(const Statistics &stat);

    inline long long getN() const   { return n; }
    inline double getMean() const   { return mean; }
    inline double getM2() const     { return M2; }

    inline double getVar() const    { return M2/static_cast<double>(n - 1); }
    inline double getCI() const     { return 1.96*std::sqrt(getVar())/std::sqrt(static_cast<double>(n)); }

protected:
    long long n;
    double mean, M2;
};
