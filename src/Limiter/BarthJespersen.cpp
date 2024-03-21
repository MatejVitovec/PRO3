#include "BarthJespersen.hpp"


double BarthJespersen::limiterFunction(double y) const
{
    return std::min(1.0, y);
}