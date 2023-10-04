#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

#include "Iawps95.hpp"

double phi0(double delta, double tau)
{
    auto val2 = [](double delta, double tau) -> double
    {
        return delta + tau;
    };

    return val2(delta, tau);
}

double T(double p0, double rho0)
{
    return 0.0;
}

int main(int argc, char** argv)
{
    double d = 10;
    double t = 5;

    Iawps95 termo = Iawps95();


    int a = 5;
}