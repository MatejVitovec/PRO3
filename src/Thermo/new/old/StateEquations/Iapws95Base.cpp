#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95Base.hpp"


Iapws95Base::Iapws95Base()
{
    loadCoeffs("Thermo/StateEquations/iapws95Coeffs/");
}

std::vector<double> Iapws95Base::loadCoeffFile(std::string name, std::string dirName, int size) const
{
    std::vector<double> out(size);
    std::ifstream f(dirName + name);
    std::string line;

    for (int i = 0; i < size; i++)
    {
        std::getline(f, line);
        out[i] = std::stod(line);
    }
            
    return out;
}

void Iapws95Base::loadCoeffs(std::string dirPath)
{
    std::vector<double> n0 = loadCoeffFile("n0", dirPath, 8);
    std::move(n0.begin(), n0.begin() + 8, coeffs0.n0.begin());
    std::vector<double> gamma0 = loadCoeffFile("gamma0", dirPath, 8);
    std::move(gamma0.begin(), gamma0.begin() + 8, coeffs0.gamma0.begin());
}


//IAPWS95BASE function

double Iapws95Base::p(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return pFunc(delta, tau);
}

double Iapws95Base::e(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return eFunc(delta, tau);
}

double Iapws95Base::s(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return sFunc(delta, tau);
}

double Iapws95Base::h(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*(critT/tau)*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Iapws95Base::a2(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    double phirdAux = phird(delta, tau);

    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Iapws95Base::a(double rho, double T) const
{
    return std::sqrt(a2(rho, T));
}

double Iapws95Base::tFromRhoE(double rho, double e, double guessT) const
{
    double delta = rho/critRho;

    double tau = nonLinearSolver.solve([=](double val) { return eFunc(delta, val) - e; },
                                       [=](double val) { return eDTauFunc(delta, val);},
                                       critT/guessT);
    
    return critT/tau;
}

double Iapws95Base::tFromRhoP(double rho, double p, double guessT) const
{
    double delta = rho/critRho;

    double tau = nonLinearSolver.solve([=](double val) { return pFunc(delta, val) - p; },
                                       [=](double val) { return pDTauFunc(delta, val);},
                                       critT/guessT);
    
    return critT/tau;
}

double Iapws95Base::rhoFromTP(double T, double p, double guessRho) const
{
    double tau = critT/T;

    double delta = nonLinearSolver.solve([=](double val) { return pFunc(val, tau) - p; },
                                         [=](double val) { return pDDeltaFunc(val, tau);},
                                         guessRho/critRho);

    return delta*critRho;
}

std::pair<double, double> Iapws95Base::rRhoFromSP(double s, double p, double guessRho, double guessT) const
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


// dimensionless Helmholtz free energy functions

double Iapws95Base::phi0(double delta, double tau) const
{
    double out = log(delta) + coeffs0.n0[0] + coeffs0.n0[1]*tau + coeffs0.n0[2]*log(tau);
    for (int i = 3; i < 8; i++)
    {
        out += coeffs0.n0[i]*log(1.0 - exp(-coeffs0.gamma0[i]*tau));
    }
    return out;
}

double Iapws95Base::phi0d(double delta, double tau) const
{
    return 1.0/delta;
}

double Iapws95Base::phi0dd(double delta, double tau) const
{
    return -1.0/(delta*delta);
}

double Iapws95Base::phi0t(double delta, double tau) const
{
    double out = coeffs0.n0[1] + coeffs0.n0[2]/tau;
    for (int i = 3; i < 8; i++)
    {
        out += coeffs0.n0[i]*coeffs0.gamma0[i]*(1/(1.0 - exp(-coeffs0.gamma0[i]*tau)) - 1.0);
    }    
    return out;
}

double Iapws95Base::phi0tt(double delta, double tau) const
{
    double out = -coeffs0.n0[2]/(tau*tau);
    for (int i = 3; i < 8; i++)
    {
        double tmp = exp(-coeffs0.gamma0[i]*tau);
        out -= coeffs0.n0[i]*coeffs0.gamma0[i]*coeffs0.gamma0[i]*tmp*pow(1.0 - tmp, -2.0);
    }    
    return out;
}

double Iapws95Base::phi0dt(double delta, double tau) const
{
    return 0.0;
}
