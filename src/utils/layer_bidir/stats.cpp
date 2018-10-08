#include "stats.h"


Statistics::Statistics(const Statistics &stats)
{
    n = stats.n;
    mean = stats.mean;
    M2 = stats.M2;
}


void Statistics::push(double x)
{
    double delta = x - mean;
    mean += delta/static_cast<double>(++n);
    M2 += delta*(x - mean);
}


void Statistics::push(const Statistics &stat)
{
    double delta = mean - stat.mean;
    mean = (mean*n + stat.mean*stat.n)/(n + stat.n);
    M2 += stat.M2 + delta*delta*n*stat.n/(n + stat.n);
    n += stat.n;
}
