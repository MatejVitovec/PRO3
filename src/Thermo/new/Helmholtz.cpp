#include <cmath>
#include <iostream>
#include <fstream>

#include "Helmholtz.hpp"

Helmholtz::Helmholtz()
{
    nonLinearSolver = NewtonMethod();
}

double Helmholtz::p(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return pFunc(delta, tau);
}

double Helmholtz::e(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return eFunc(delta, tau);
}

double Helmholtz::s(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return sFunc(delta, tau);
}

double Helmholtz::h(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*(critT/tau)*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Helmholtz::a2(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    double phirdAux = phird(delta, tau);

    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Helmholtz::a(double rho, double T) const
{
    return std::sqrt(a2(rho, T));
}

double Helmholtz::tFromRhoE(double rho, double e, double guessT) const
{
    double delta = rho/critRho;

    double tau = nonLinearSolver.solve([=](double val) { return eFunc(delta, val) - e; },
                                       [=](double val) { return eDTauFunc(delta, val);},
                                       critT/guessT);
    
    return critT/tau;
}

double Helmholtz::tFromRhoP(double rho, double p, double guessT) const
{
    double delta = rho/critRho;

    double tau = nonLinearSolver.solve([=](double val) { return pFunc(delta, val) - p; },
                                       [=](double val) { return pDTauFunc(delta, val); },
                                       critT/guessT);
    
    return critT/tau;
}

double Helmholtz::rhoFromTP(double T, double p, double guessRho) const
{
    double tau = critT/T;

    double delta = nonLinearSolver.solve([=](double val) { return pFunc(val, tau) - p; },
                                         [=](double val) { return pDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}

std::pair<double, double> Helmholtz::rRhoFromSP(double s, double p, double guessRho, double guessT) const
{
    std::pair<double, double> result = nonLinearSolver.solve([=](double val1, double val2) { return sFunc(val1, val2) - s; },
                                                             [=](double val1, double val2) { return pFunc(val1, val2) - p; },
                                                             sDDeltaFunc,
                                                             sDTauFunc,
                                                             pDDeltaFunc,
                                                             pDTauFunc,
                                                             guessRho/critRho,
                                                             critT/guessT);

    return std::make_pair(result.first*critRho, critT/result.second);
}
