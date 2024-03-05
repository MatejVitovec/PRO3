#include "CubicLimiter.hpp"

double CubicLimiter::limiterFunction(double y) const
{
    if(y < K)
    {
        return a*std::pow(y, 3) + b*std::pow(y, 2) + y;
    }
    else
    {
        return 1.0;
    }
        
}