#include <cmath>
#include <iostream>
#include <fstream>

#include "SpecialGas.hpp"


SpecialGas::SpecialGas() : Helmholtz()
{
    loadCoeffs("Thermo/StateEquations/specialGasCoeffs/");
}

std::vector<double> SpecialGas::loadCoeffFile(std::string name, std::string dirName, int size) const
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

void SpecialGas::loadCoeffs(std::string dirPath)
{
    std::vector<double> n0 = loadCoeffFile("n0", dirPath, 8);
    std::move(n0.begin(), n0.begin() + 8, coeffs.n0.begin());
    std::vector<double> gamma0 = loadCoeffFile("gamma0", dirPath, 8);
    std::move(gamma0.begin(), gamma0.begin() + 8, coeffs.gamma0.begin());    
    std::vector<double> d = loadCoeffFile("d", dirPath, 7);
    std::move(d.begin(), d.begin() + 7, coeffs.d.begin());
    std::vector<double> t = loadCoeffFile("t", dirPath, 7);
    std::move(t.begin(), t.begin() + 7, coeffs.t.begin());
    std::vector<double> n = loadCoeffFile("n", dirPath, 7);
    std::move(n.begin(), n.begin() + 7, coeffs.n.begin());
}

// dimensionless Helmholtz free energy functions

double SpecialGas::phi0(double delta, double tau) const
{
    double out = log(delta) + coeffs.n0[0] + coeffs.n0[1]*tau + coeffs.n0[2]*log(tau);
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*log(1.0 - exp(-coeffs.gamma0[i]*tau));
    }
    return out;
}

double SpecialGas::phi0d(double delta, double tau) const
{
    return 1.0/delta;
}

double SpecialGas::phi0dd(double delta, double tau) const
{
    return -1.0/(delta*delta);
}

double SpecialGas::phi0t(double delta, double tau) const
{
    double out = coeffs.n0[1] + coeffs.n0[2]/tau;
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*coeffs.gamma0[i]*(1/(1.0 - exp(-coeffs.gamma0[i]*tau)) - 1.0);
    }    
    return out;
}

double SpecialGas::phi0tt(double delta, double tau) const
{
    double out = -coeffs.n0[2]/(tau*tau);
    for (int i = 3; i < 8; i++)
    {
        double tmp = exp(-coeffs.gamma0[i]*tau);
        out -= coeffs.n0[i]*coeffs.gamma0[i]*coeffs.gamma0[i]*tmp*pow(1.0 - tmp, -2.0);
    }    
    return out;
}

double SpecialGas::phi0dt(double delta, double tau) const
{
    return 0.0;
}

double SpecialGas::phir(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i]);
    }
    return out;
}

double SpecialGas::phird(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i]);
    }
    return out;
}

double SpecialGas::phirdd(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*(coeffs.d[i] - 1.0)*pow(delta, coeffs.d[i] - 2.0)*pow(tau, coeffs.t[i]);
    }
    return out;
}

double SpecialGas::phirt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 1.0);
    }
    return out;
}

double SpecialGas::phirtt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*(coeffs.t[i] - 1.0)*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 2.0);
    }
    return out;
}

double SpecialGas::phirdt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*coeffs.t[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i] - 1.0);
    }
    return out;
}