#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

#include "Iawps95.hpp"


int main(int argc, char** argv)
{
    Iawps95 termo = Iawps95();

    double T = 500.0;
    double rho = 1084.564;
    //double e = 1966.95;

    //double T = termo.implicitTemperature(rho, e, 600.0);
    
    double p = termo.p(rho, T);
    double e = termo.e(rho, T);
    double s = termo.s(rho, T);
    double w = sqrt(termo.w2(rho, T));

    //termo.test(838.025, 500);



    int a = 5;
}