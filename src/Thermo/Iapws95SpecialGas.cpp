#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95SpecialGas.hpp"


Iapws95SpecialGas::Iapws95SpecialGas(): Thermo()
{
    loadCoeffs("Thermo/iapws95SpecialGasCoeffs/");
}

std::vector<double> Iapws95SpecialGas::loadCoeffFile(std::string name, std::string dirName, int size) const
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

void Iapws95SpecialGas::loadCoeffs(std::string dirPath)
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

Vars<3> Iapws95SpecialGas::updateThermo(const Compressible& data, const Compressible& dataOld) const
{
    //treti clen stara hodnota teploty jako odhad pro nelinearni reseni rovnice
    double rho = data.density();
    double T = tFromRhoE(rho, data.internalEnergy(), dataOld.temperature());

    return Vars<3>({T, p(rho, T), a(rho, T)});
}

Compressible Iapws95SpecialGas::primitiveToConservative(const Vars<5>& primitive) const
{
    double rho = primitive[0];
    double p = primitive[4];
    //odhad T pomoci idealniho plynu
    double T = tFromRhoP(rho, p, p/(specGasConst*rho));

    return Compressible({rho,
                         rho*primitive[1],
                         rho*primitive[2],
                         rho*primitive[3],
                         rho*(e(rho, T) + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {T, p, a(rho, T)});
}

Compressible Iapws95SpecialGas::isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    double pIn = std::min(stateIn.pressure(), pTot);

    double sTot = s(rhoTot, TTot);

    //newton method for 2 equations
    double guessRho = stateIn.density();
    double guessT = stateIn.temperature();

    double deltaOld = guessRho/critRho;
    double tauOld = critT/guessT;
    double delta = 0.0;
    double tau = 0.0;

    double change = 10000000.0;
    int i = 0;

    // Newton method loop
    while (change > numericalTolerance)
    {
        double auxPhird = phird(deltaOld, tauOld);
        double auxPhirdt = phirdt(deltaOld, tauOld);
        double auxConst = specGasConst*critRho*critT*(deltaOld/tauOld);

        double f = sTot - specGasConst*(tauOld*(phi0t(deltaOld, tauOld) + phirt(deltaOld, tauOld)) - phi0(deltaOld, tauOld) - phir(deltaOld, tauOld));
        double g = pIn - auxConst - auxConst*deltaOld*auxPhird;
        double fd = -specGasConst*(tauOld*(phi0dt(deltaOld, tauOld) + auxPhirdt) - phi0d(deltaOld, tauOld) - auxPhird);
        double ft = -specGasConst*tauOld*(phi0tt(deltaOld, tauOld) + phirtt(deltaOld, tauOld));
        double gd = -specGasConst*critRho*(critT/tauOld) - 2.0*auxConst*auxPhird - auxConst*deltaOld*phirdd(deltaOld, tauOld);
        double gt = auxConst/tauOld + auxConst*(deltaOld/tauOld)*auxPhird - auxConst*deltaOld*auxPhirdt;
        
        delta = deltaOld - (f*gt - g*ft)/(fd*gt - ft*gd);
        tau = tauOld - (g*fd - f*gd)/(fd*gt - ft*gd);
    
        change = std::max(fabs((delta - deltaOld)/delta), fabs((tau - tauOld)/tau));

        deltaOld = delta;
        tauOld = tau;
        i++;
    }

    double rho = delta*critRho;
    double T = critT/tau;

    double absU2 = std::max(2.0*h(rhoTot, TTot) - 2.0*h(rho, T), 0.0);
    double absU = std::sqrt(absU2);


    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e(rho, T))},
                         {T, pIn, a(rho, T)});
}

Compressible Iapws95SpecialGas::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

Compressible Iapws95SpecialGas::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double TTot = tFromRhoP(rhoTot, pTot, pTot/(specGasConst*rhoTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

double Iapws95SpecialGas::p(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return rho*specGasConst*T*(1.0 + delta*phird(delta, tau));
}

double Iapws95SpecialGas::e(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau));
}

double Iapws95SpecialGas::s(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau));
}

double Iapws95SpecialGas::h(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*T*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Iapws95SpecialGas::a2(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    double phirdAux = phird(delta, tau);

    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Iapws95SpecialGas::a(double rho, double T) const
{
    return std::sqrt(a2(rho, T));
}

double Iapws95SpecialGas::tFromRhoE(double rho, double e, double guessT) const
{
    double delta = rho/critRho;
    double tauOld = critT/guessT;
    double tau = 0.0;

    double change = 100000;
    int i = 0;

    // Newton method
    while (change > numericalTolerance)
    {
        tau = tauOld - (e - specGasConst*critT*(phi0t(delta, tauOld) + phirt(delta, tauOld)))/(-specGasConst*critT*(phi0tt(delta, tauOld) + phirtt(delta, tauOld)));
    
        change = fabs((tau - tauOld)/*/tau*/);

        tauOld = tau;
        i++;
    }
    
    return critT/tau;
}

double Iapws95SpecialGas::tFromRhoP(double rho, double p, double guessT) const
{
    double delta = rho/critRho;
    double tauOld = critT/guessT;
    double tau = 0.0;

    double change = 100000;
    int i = 0;

    // Newton method
    while (change > numericalTolerance)
    {
        double phirdAux = phird(delta, tauOld);
        double auxVar = specGasConst*critRho*critT*(delta/tauOld);
        tau = tauOld - (p - auxVar*(1.0 + delta*phirdAux))/(auxVar*(1.0/tauOld + (delta/tauOld)*phirdAux - delta*phirdt(delta, tauOld)));
        
        change = fabs((tau - tauOld)/*/tau*/);

        tauOld = tau;
        i++;
    }
    
    return critT/tau;
}

double Iapws95SpecialGas::rhoFromTP(double T, double p, double guessRho) const
{
    double tau = critT/T;
    double deltaOld = guessRho/critRho;    
    double delta = 0.0;

    double change = 100000;
    int i = 0;

    // Newton method
    while (change > numericalTolerance)
    {
        double auxVar = critRho*critT*(specGasConst/tau);
        double phirdAux = phird(deltaOld, tau);
        delta = deltaOld - (p - auxVar*(deltaOld + deltaOld*deltaOld*phirdAux))/(auxVar*(1.0 + 2.0*deltaOld*phirdAux + deltaOld*deltaOld*phirdt(deltaOld, tau)));
        
        change = fabs((delta - deltaOld)/*/delta*/);

        deltaOld = delta;
        i++;
    }
    
    return delta*critRho;
}



// dimensionless Helmholtz free energy functions

double Iapws95SpecialGas::phi0(double delta, double tau) const
{
    double out = log(delta) + coeffs.n0[0] + coeffs.n0[1]*tau + coeffs.n0[2]*log(tau);
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*log(1.0 - exp(-coeffs.gamma0[i]*tau));
    }
    return out;
}

double Iapws95SpecialGas::phi0d(double delta, double tau) const
{
    return 1.0/delta;
}

double Iapws95SpecialGas::phi0dd(double delta, double tau) const
{
    return -1.0/(delta*delta);
}

double Iapws95SpecialGas::phi0t(double delta, double tau) const
{
    double out = coeffs.n0[1] + coeffs.n0[2]/tau;
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*coeffs.gamma0[i]*(1/(1.0 - exp(-coeffs.gamma0[i]*tau)) - 1.0);
    }    
    return out;
}

double Iapws95SpecialGas::phi0tt(double delta, double tau) const
{
    double out = -coeffs.n0[2]/(tau*tau);
    for (int i = 3; i < 8; i++)
    {
        double tmp = exp(-coeffs.gamma0[i]*tau);
        out -= coeffs.n0[i]*coeffs.gamma0[i]*coeffs.gamma0[i]*tmp*pow(1.0 - tmp, -2.0);
    }    
    return out;
}

double Iapws95SpecialGas::phi0dt(double delta, double tau) const
{
    return 0.0;
}

double Iapws95SpecialGas::phir(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i]);
    }
    return out;
}

double Iapws95SpecialGas::phird(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i]);
    }
    return out;
}

double Iapws95SpecialGas::phirdd(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*(coeffs.d[i] - 1.0)*pow(delta, coeffs.d[i] - 2.0)*pow(tau, coeffs.t[i]);
    }
    return out;
}

double Iapws95SpecialGas::phirt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 1.0);
    }
    return out;
}

double Iapws95SpecialGas::phirtt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*(coeffs.t[i] - 1.0)*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 2.0);
    }
    return out;
}

double Iapws95SpecialGas::phirdt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*coeffs.t[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i] - 1.0);
    }
    return out;
}
