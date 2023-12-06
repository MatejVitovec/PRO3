#include <cmath>

#include "NewtonMethod.hpp"

double NewtonMethod::solve(std::function<double(double)> f, std::function<double(double)> df, double guess) const
{
    double val = guess;
    double valNew = 0.0;

    double change = 100000;
    int i = 0;

    while (change > numTol)
    {
        valNew = val - f(val)/df(val);
    
        change = fabs((valNew - val)/valNew);

        val = valNew;
        i++;
    }
    
    return val;
}

std::pair<double, double> NewtonMethod::solve(std::function<double(double, double)> f,
                                              std::function<double(double, double)> fd1,
                                              std::function<double(double, double)> fd2,
                                              std::function<double(double, double)> g,
                                              std::function<double(double, double)> gd1,
                                              std::function<double(double, double)> gd2,
                                              double guess1,
                                              double guess2) const
{
    double val1 = guess1;
    double val2 = guess2;
    double val1New = 0.0;
    double val2New = 0.0;

    double change = 10000000.0;
    int i = 0;

    // Newton method loop
    while (change > numTol)
    {        
        val1New = val1 - (f(val1, val2)*gd2(val1, val2) - g(val1, val2)*fd2(val1, val2))/(fd1(val1, val2)*gd2(val1, val2) - fd2(val1, val2)*gd1(val1, val2));
        val2New = val2 - (g(val1, val2)*fd1(val1, val2) - f(val1, val2)*gd1(val1, val2))/(fd1(val1, val2)*gd2(val1, val2) - fd2(val1, val2)*gd1(val1, val2));
    
        change = std::max(fabs((val1New - val1)/val1New), fabs((val2New - val2)/val2New));

        val1 = val1New;
        val2 = val2New;
        i++;
    }

    return std::make_pair(val1, val2);
}