#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95.hpp"


Iapws95::Iapws95(): Thermo()
{
    loadCoeffs("Thermo/iapws95Coeffs/");
}

std::vector<double> Iapws95::loadCoeffFile(std::string name, std::string dirName, int size) const
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

void Iapws95::loadCoeffs(std::string dirPath)
{
    std::vector<double> n0 = loadCoeffFile("n0", dirPath, 8);
    std::move(n0.begin(), n0.begin() + 8, coeffs.n0.begin());
    std::vector<double> gamma0 = loadCoeffFile("gamma0", dirPath, 8);
    std::move(gamma0.begin(), gamma0.begin() + 8, coeffs.gamma0.begin());
    std::vector<double> c = loadCoeffFile("c", dirPath, 56);
    std::move(c.begin(), c.begin() + 56, coeffs.c.begin());
    std::vector<double> d = loadCoeffFile("d", dirPath, 56);
    std::move(d.begin(), d.begin() + 56, coeffs.d.begin());
    std::vector<double> t = loadCoeffFile("t", dirPath, 56);
    std::move(t.begin(), t.begin() + 56, coeffs.t.begin());
    std::vector<double> n = loadCoeffFile("n", dirPath, 56);
    std::move(n.begin(), n.begin() + 56, coeffs.n.begin());
    std::vector<double> alpha = loadCoeffFile("alpha", dirPath, 56);
    std::move(alpha.begin(), alpha.begin() + 56, coeffs.alpha.begin());
    std::vector<double> beta = loadCoeffFile("beta", dirPath, 56);
    std::move(beta.begin(), beta.begin() + 56, coeffs.beta.begin());
    std::vector<double> gamma = loadCoeffFile("gamma", dirPath, 56);
    std::move(gamma.begin(), gamma.begin() + 56, coeffs.gamma.begin());
    std::vector<double> epsilon = loadCoeffFile("epsilon", dirPath, 56);
    std::move(epsilon.begin(), epsilon.begin() + 56, coeffs.epsilon.begin());
    std::vector<double> a = loadCoeffFile("a", dirPath, 56);
    std::move(a.begin(), a.begin() + 56, coeffs.a.begin());
    std::vector<double> b = loadCoeffFile("b", dirPath, 56);
    std::move(b.begin(), b.begin() + 56, coeffs.b.begin());
    std::vector<double> B = loadCoeffFile("B", dirPath, 56);
    std::move(B.begin(), B.begin() + 56, coeffs.B.begin());
    std::vector<double> C = loadCoeffFile("C", dirPath, 56);
    std::move(C.begin(), C.begin() + 56, coeffs.C.begin());
    std::vector<double> D = loadCoeffFile("D", dirPath, 56);
    std::move(D.begin(), D.begin() + 56, coeffs.D.begin());
    std::vector<double> A = loadCoeffFile("A", dirPath, 56);
    std::move(A.begin(), A.begin() + 56, coeffs.A.begin());
}

Vars<3> Iapws95::updateThermo(const Compressible& data, const Compressible& dataOld) const
{
    //treti clen stara hodnota teploty jako odhad pro nelinearni reseni rovnice
    double rho = data.density();
    double T = tFromRhoE(rho, data.internalEnergy(), dataOld.temperature());

    return Vars<3>({T, p(rho, T), a(rho, T)});
}

Compressible Iapws95::primitiveToConservative(const Vars<5>& primitive) const
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

Compressible Iapws95::isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
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

Compressible Iapws95::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

Compressible Iapws95::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double TTot = tFromRhoP(rhoTot, pTot, pTot/(specGasConst*rhoTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

double Iapws95::p(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return rho*specGasConst*T*(1.0 + delta*phird(delta, tau));
}

double Iapws95::e(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau));
}

double Iapws95::s(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau));
}

double Iapws95::h(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*T*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Iapws95::a2(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    double phirdAux = phird(delta, tau);

    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Iapws95::a(double rho, double T) const
{
    return std::sqrt(a2(rho, T));
}

double Iapws95::tFromRhoE(double rho, double e, double guessT) const
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

double Iapws95::tFromRhoP(double rho, double p, double guessT) const
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

/*double Iapws95::tFromRhoS(double rho, double s, double guessT) const
{
    double delta = rho/critRho;
    double tauOld = critT/guessT;
    double tau = 0.0;

    double change = 100000;
    int i = 0;

    // Newton method
    while (change > numericalTolerance)
    {
        tau = tauOld - (s/specGasConst - tauOld*(phi0t(delta, tauOld) + phirt(delta, tauOld)) + phi0(delta, tauOld) + phir(delta, tauOld))/(-tauOld*(phi0tt(delta, tauOld) + phirtt(delta, tauOld)));
        
        change = fabs((tau - tauOld)/tau);

        tauOld = tau;
        i++;
    }
    
    return critT/tau;
}*/

double Iapws95::rhoFromTP(double T, double p, double guessRho) const
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

double Iapws95::phi0(double delta, double tau) const
{
    double out = log(delta) + coeffs.n0[0] + coeffs.n0[1]*tau + coeffs.n0[2]*log(tau);
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*log(1.0 - exp(-coeffs.gamma0[i]*tau));
    }
    return out;
}

double Iapws95::phi0d(double delta, double tau) const
{
    return 1.0/delta;
}

double Iapws95::phi0dd(double delta, double tau) const
{
    return -1.0/(delta*delta);
}

double Iapws95::phi0t(double delta, double tau) const
{
    double out = coeffs.n0[1] + coeffs.n0[2]/tau;
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*coeffs.gamma0[i]*(1/(1.0 - exp(-coeffs.gamma0[i]*tau)) - 1.0);
    }    
    return out;
}

double Iapws95::phi0tt(double delta, double tau) const
{
    double out = -coeffs.n0[2]/(tau*tau);
    for (int i = 3; i < 8; i++)
    {
        double tmp = exp(-coeffs.gamma0[i]*tau);
        out -= coeffs.n0[i]*coeffs.gamma0[i]*coeffs.gamma0[i]*tmp*pow(1.0 - tmp, -2.0);
    }    
    return out;
}

double Iapws95::phi0dt(double delta, double tau) const
{
    return 0.0;
}

double Iapws95::phir(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i]);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-pow(delta, coeffs.c[i]));
    }
    for (int i = 51; i < 54; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*pow(delta - coeffs.epsilon[i], 2) - coeffs.beta[i]*pow(tau - coeffs.gamma[i], 2));
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i])*delta*psiFunc(delta, tau, i);
    }
    return out;
}

double Iapws95::phird(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i]);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*exp(-pow(delta, coeffs.c[i]))*(pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i])*(coeffs.d[i] - coeffs.c[i]*pow(delta, coeffs.c[i])));
    }
    for (int i = 51; i < 54; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*pow(delta - coeffs.epsilon[i], 2) - coeffs.beta[i]*pow(tau - coeffs.gamma[i], 2));
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*(pow(deltaFunc(delta, tau, i), coeffs.b[i])*(psiFunc(delta, tau, i) + delta*psiFuncd(delta, tau, i)) + deltaFuncbid(delta, tau, i)*delta*psiFunc(delta, tau, i));                            
    }
    return out;
}

double Iapws95::phirdd(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*(coeffs.d[i] - 1.0)*pow(delta, coeffs.d[i] - 2.0)*pow(tau, coeffs.t[i]);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*exp(-pow(delta, coeffs.c[i]))*(pow(delta, coeffs.d[i] - 2.0)*pow(tau, coeffs.t[i])*((coeffs.d[i] - coeffs.c[i]*pow(delta, coeffs.c[i]))*(coeffs.d[i] - 1.0 - coeffs.c[i]*pow(delta, coeffs.c[i])) - pow(coeffs.c[i], 2)*pow(delta, coeffs.c[i])));
    }
    for (int i = 51; i < 54; i++)
    {
        out += coeffs.n[i]*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*pow(delta - coeffs.epsilon[i], 2) - coeffs.beta[i]*pow(tau - coeffs.gamma[i], 2))*(-2.0*coeffs.alpha[i]*pow(delta, coeffs.d[i]) + 4.0*pow(coeffs.alpha[i], 2)*pow(delta, coeffs.d[i])*pow(delta - coeffs.epsilon[i], 2) - 4.0*coeffs.d[i]*coeffs.alpha[i]*pow(delta, coeffs.d[i] - 1.0)*(delta - coeffs.epsilon[i]) + coeffs.d[i]*(coeffs.d[i] - 1.0)*pow(delta, coeffs.d[i] - 2.0));
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*(pow(deltaFunc(delta, tau, i), coeffs.b[i])*(2*psiFuncd(delta, tau, i) + delta*psiFuncdd(delta, tau, i)) + 2*deltaFuncbid(delta, tau, i)*(psiFunc(delta, tau, i) + delta*psiFunc(delta, tau ,i)) + deltaFuncbidd(delta, tau, i)*delta*psiFunc(delta, tau, i));
    }
    return out;
}

double Iapws95::phirt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 1.0);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 1.0)*exp(-pow(delta, coeffs.c[i]));
    }
    for (int i = 51; i < 54; i++)
    {
        out += (coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*(delta - coeffs.epsilon[i])*(delta - coeffs.epsilon[i]) - coeffs.beta[i]*(tau - coeffs.gamma[i])*(tau - coeffs.gamma[i])))*(coeffs.t[i]/tau - 2.0*coeffs.beta[i]*(tau - coeffs.gamma[i]));                           
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*delta*(deltaFuncbit(delta, tau, i)*psiFunc(delta, tau, i) + pow(deltaFunc(delta, tau, i), coeffs.b[i])*psiFunct(delta, tau, i));
    }
    return out;
}

double Iapws95::phirtt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*(coeffs.t[i] - 1.0)*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 2.0);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*(coeffs.t[i] - 1.0)*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i] - 2.0)*exp(-pow(delta, coeffs.c[i]));
    }
    for (int i = 51; i < 54; i++)
    {
        out += (coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*(delta - coeffs.epsilon[i])*(delta - coeffs.epsilon[i]) - coeffs.beta[i]*(tau - coeffs.gamma[i])*(tau - coeffs.gamma[i])))*(pow(coeffs.t[i]/tau - 2.0*coeffs.beta[i]*(tau - coeffs.gamma[i]), 2.0) - coeffs.t[i]/(tau*tau) - 2.0*coeffs.beta[i]);                           
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*delta*(deltaFuncbitt(delta, tau, i)*psiFunc(delta, tau, i) + 2.0*deltaFuncbidt(delta, tau, i)*psiFunct(delta, tau, i) + pow(deltaFunc(delta, tau, i), coeffs.b[i])*psiFunctt(delta, tau, i));
    }
    return out;
}

double Iapws95::phirdt(double delta, double tau) const
{
    double out = 0.0;
    for (int i = 0; i < 7; i++)
    {
        out += coeffs.n[i]*coeffs.d[i]*coeffs.t[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i] - 1.0);
    }
    for (int i = 7; i < 51; i++)
    {
        out += coeffs.n[i]*coeffs.t[i]*pow(delta, coeffs.d[i] - 1.0)*pow(tau, coeffs.t[i] - 1.0)*(coeffs.d[i] - coeffs.c[i]*pow(delta, coeffs.c[i]))*exp(-pow(delta, coeffs.c[i]));
    }
    for (int i = 51; i < 54; i++)
    {
        out += coeffs.n[i]*pow(delta, coeffs.d[i])*pow(tau, coeffs.t[i])*exp(-coeffs.alpha[i]*pow(delta - coeffs.epsilon[i], 2) - coeffs.beta[i]*pow(tau - coeffs.gamma[i], 2))*(coeffs.d[i]/delta - 2.0*coeffs.alpha[i]*(delta - coeffs.epsilon[i]))*(coeffs.t[i]/tau - 2.0*coeffs.beta[i]*(tau - coeffs.gamma[i]));                           
    }
    for (int i = 54; i < 56; i++)
    {
        out += coeffs.n[i]*(pow(deltaFunc(delta, tau, i), coeffs.b[i])*(psiFunct(delta, tau, i) + delta*psiFuncdt(delta, tau, i)) + delta*(deltaFuncbid(delta, tau, i)*psiFunct(delta, tau, i)) + deltaFuncbit(delta, tau, i)*(psiFunct(delta, tau, i) + delta*psiFuncd(delta, tau, i)) + deltaFuncbidt(delta, tau, i)*delta*psiFunc(delta, tau, i));
    }
    return out;
}


//aux functions

double Iapws95::deltaFunc(double delta, double tau, int i) const
{
    return pow(thetaFunc(delta, tau, i), 2) + coeffs.B[i]*pow(pow(delta - 1.0, 2), coeffs.a[i]);
}

double Iapws95::thetaFunc(double delta, double tau, int i) const
{
    return 1.0 - tau + coeffs.A[i]*pow(pow(delta - 1.0, 2), 1.0/(2.0*coeffs.beta[i]));
}

double Iapws95::psiFunc(double delta, double tau, int i) const
{
    return exp(-coeffs.C[i]*pow(delta - 1.0, 2) - coeffs.D[i]*pow(tau - 1.0, 2));
}


double Iapws95::deltaFuncd(double delta, double tau, int i) const
{
    return (delta - 1.0)*(coeffs.A[i]*thetaFunc(delta, tau, i)*(2.0/coeffs.beta[i])*pow(pow(delta - 1.0, 2), 1/(2.0*coeffs.beta[i]) - 1.0) + 2.0*coeffs.B[i]*coeffs.a[i]*pow(pow(delta - 1.0, 2.0), coeffs.a[i] - 1.0));
}

double Iapws95::deltaFuncdd(double delta, double tau, int i) const
{
    return (1.0/(delta - 1.0))*deltaFuncd(delta, tau, i) + pow(delta - 1.0, 2)*(4.0*coeffs.B[i]*coeffs.a[i]*(coeffs.a[i] - 1.0)*pow(pow(delta - 1.0, 2), coeffs.a[i] - 2.0) + 2.0*pow(coeffs.A[i], 2)*pow(1/coeffs.beta[i], 2)*pow(pow(pow(delta - 1.0, 2), 1/(2*coeffs.beta[i]) - 1.0), 2) + coeffs.A[i]*thetaFunc(delta, tau, i)*(4.0/coeffs.beta[i])*(1/(2*coeffs.beta[i]) - 1.0)*pow(pow(delta - 1.0, 2), 1/(2*coeffs.beta[i] - 2.0)));
}


double Iapws95::deltaFuncbid(double delta, double tau, int i) const
{
    return coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*deltaFuncd(delta, tau, i);
}

double Iapws95::deltaFuncbidd(double delta, double tau, int i) const
{
    return coeffs.b[i]*(pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*deltaFuncdd(delta, tau, i) + (coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0)*pow(deltaFuncd(delta, tau, i), 2));
}

double Iapws95::deltaFuncbit(double delta, double tau, int i) const
{
    return -2.0*thetaFunc(delta, tau, i)*coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0);
}

double Iapws95::deltaFuncbitt(double delta, double tau, int i) const
{
    return 2.0*coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0) + 4.0*pow(thetaFunc(delta, tau, i), 2.0)*coeffs.b[i]*(coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0);
}

double Iapws95::deltaFuncbidt(double delta, double tau, int i) const
{
    return -coeffs.A[i]*coeffs.b[i]*(2.0/coeffs.beta[i])*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*(delta - 1.0)*pow(pow(delta - 1, 2), 1.0/(2.0*coeffs.beta[i]) - 1.0) - 2.0*thetaFunc(delta, tau, i)*coeffs.b[i]*(coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0)*deltaFuncd(delta, tau, i);
}


double Iapws95::psiFuncd(double delta, double tau, int i) const
{
    return -2.0*coeffs.C[i]*(delta - 1.0)*psiFunc(delta, tau, i);
}

double Iapws95::psiFuncdd(double delta, double tau, int i) const
{
    return (2.0*coeffs.C[i]*pow(delta - 1.0, 2) - 1.0)*2.0*coeffs.C[i]*psiFunc(delta, tau, i);
}

double Iapws95::psiFunct(double delta, double tau, int i) const
{
    return -2.0*coeffs.D[i]*(tau - 1.0)*psiFunc(delta, tau, i);
}

double Iapws95::psiFunctt(double delta, double tau, int i) const
{
    return (2.0*coeffs.D[i]*(tau - 1.0)*(tau - 1.0) - 1.0)*2.0*coeffs.D[i]*psiFunc(delta, tau, i);
}

double Iapws95::psiFuncdt(double delta, double tau, int i) const
{
    return 4.0*coeffs.C[i]*coeffs.D[i]*(delta - 1.0)*(tau - 1.0)*psiFunc(delta, tau, i);
}
