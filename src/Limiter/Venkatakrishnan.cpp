#include "Venkatakrishnan.hpp"
#include <iostream>


double Venkatakrishnan::limiterFunction(double y) const
{
    return (y*y + 2*y)/(y*y + y + 2);
}