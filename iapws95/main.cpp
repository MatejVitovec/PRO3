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

    //double T = 500.0;
    //double rho = 1084.564;

    //double rho = 1084.564;
    //double s = 2032.3750919;
    //double T = 500.0;

    //double rho = 1080;
    //double e = 768.08;
    
    //double p = 98071.7;
    //double T = 298.65;

    //double T = termo.temperatureFromRhoE(rho, e, 300.0);
    //double T = termo.temperatureFromRhoP(rho, p, 300.0);

    //double TIG = e*(1.4 - 1.0)/0.46151805;
    
    /*double p = termo.p(rho, T);
    double e = termo.e(rho, T);
    double s = termo.s(rho, T);
    double w = sqrt(termo.w2(rho, T));*/

    //double rho = termo.implicitRhoFromTS(T, s, 800.0);

    //termo.test(838.025, 500);

    termo.saturatedRhoEData();


    int a = 5;
}