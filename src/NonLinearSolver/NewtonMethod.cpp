#include <cmath>

#include "NewtonMethod.hpp"

double NewtonMethod::solve(std::function<double(double)> f, std::function<double(double)> df, double guess) const
{
    double valOld = guess;
    double val = 0.0;

    double change = 100000;
    int i = 0;

    while (change > numTol)
    {
        val = valOld - f(val)/df(val);
    
        change = fabs((val - valOld)/val);

        valOld = val;
        i++;
    }
    
    return val;
}

std::pair<double, double> NewtonMethod::solve(std::function<double(double, double)> f, std::function<double(double, double)> fd1, std::function<double(double, double)> fd2,
                                               std::function<double(double, double)> g, std::function<double(double, double)> gd1, std::function<double(double, double)> gd2,
                                               double guess1, double guess2) const
{
    double val1Old = guess1;
    double val2Old = guess2;
    double val1 = 0.0;
    double val2 = 0.0;

    double change = 10000000.0;
    int i = 0;

    // Newton method loop
    while (change > numTol)
    {        
        val1 = val1Old - (f(val1, val2)*gd2(val1, val2) - g(val1, val2)*fd2(val1, val2))/(fd1(val1, val2)*gd2(val1, val2) - fd2(val1, val2)*gd1(val1, val2));
        val2 = val2Old - (g(val1, val2)*fd1(val1, val2) - f(val1, val2)*gd1(val1, val2))/(fd1(val1, val2)*gd2(val1, val2) - fd2(val1, val2)*gd1(val1, val2));
    
        change = std::max(fabs((val1 - val1Old)/val1), fabs((val2 - val2Old)/val2));

        val1Old = val1;
        val2Old = val2;
        i++;
    }

    return std::make_pair(val1, val2);
}