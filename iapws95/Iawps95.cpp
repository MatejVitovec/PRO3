#include <cmath>

#include "Iawps95.hpp"


Iawps95::Iawps95()
{
    loadCoeffs("coeffs/");
}

std::vector<double> Iawps95::loadCoeffFile(std::string name, std::string dirName, int size) const
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

void Iawps95::loadCoeffs(std::string dirPath)
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

void Iawps95::test(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    std::cout << "phi0: " << phi0(delta, tau) << std::endl;
    std::cout << "phi0d: " << phi0d(delta, tau) << std::endl;
    std::cout << "phi0dd: " << phi0dd(delta, tau) << std::endl;
    std::cout << "phi0t: " << phi0t(delta, tau) << std::endl;
    std::cout << "phi0tt: " << phi0tt(delta, tau) << std::endl;
    std::cout << "phi0dt: " << phi0dt(delta, tau) << std::endl;

    std::cout << std::endl;

    std::cout << "phir: " << phir(delta, tau) << std::endl;
    std::cout << "phird: " << phird(delta, tau) << std::endl;
    std::cout << "phirdd: " << phirdd(delta, tau) << std::endl;
    std::cout << "phirt: " << phirt(delta, tau) << std::endl;
    std::cout << "phirtt: " << phirtt(delta, tau) << std::endl;
    std::cout << "phirdt: " << phirdt(delta, tau) << std::endl;
}



double Iawps95::p(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return rho*specGasConst*T*(1.0 + delta*phird(delta, tau));
}

double Iawps95::e(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau));
}

double Iawps95::s(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau));
}

double Iawps95::h(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    return specGasConst*T*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Iawps95::w2(double rho, double T) const
{
    double delta = rho/critRho;
    double tau = critT/T;

    double phirdAux = phird(delta, tau);

    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Iawps95::implicitTFromRhoE(double rho, double e, double guessT) const
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
    
        change = fabs(tau - tauOld);

        tauOld = tau;
        i++;
    }

    //std::cout << i << std::endl;
    
    return critT/tau;
}

double Iawps95::implicitTFromRhoP(double rho, double p, double guessT) const
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
        
        change = fabs(tau - tauOld);

        tauOld = tau;
        i++;
    }

    std::cout << i << std::endl;
    
    return critT/tau;
}

double Iawps95::saturatedVaporDensity(double T) const
{
    double stauratedCoeffs[6] = {-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063};
    double theta = 1.0 - T/critT;

    return critRho*exp(stauratedCoeffs[0]*pow(theta, 2.0/6.0) + 
                       stauratedCoeffs[1]*pow(theta, 4.0/6.0) + 
                       stauratedCoeffs[2]*pow(theta, 8.0/6.0) + 
                       stauratedCoeffs[3]*pow(theta, 18.0/6.0) + 
                       stauratedCoeffs[4]*pow(theta, 37.0/6.0) + 
                       stauratedCoeffs[5]*pow(theta, 71.0/6.0));
}

double Iawps95::saturatedLiquidDensity(double T) const
{
    double stauratedCoeffs[6] = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.450};
    double theta = 1.0 - T/critT;

    return critRho*(1.0 +
                    stauratedCoeffs[0]*pow(theta, 1.0/3.0) + 
                    stauratedCoeffs[1]*pow(theta, 2.0/3.0) + 
                    stauratedCoeffs[2]*pow(theta, 5.0/3.0) + 
                    stauratedCoeffs[3]*pow(theta, 16.0/3.0) + 
                    stauratedCoeffs[4]*pow(theta, 43.0/3.0) + 
                    stauratedCoeffs[5]*pow(theta, 110.0/3.0));
}

double Iawps95::vaporPressure(double T) const
{
    double vaporCoeffs[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
    double critP = 22064000.0;
    double theta = 1.0 - T/critT;

    return critP*exp(critT/T*(vaporCoeffs[0]*theta + 
                              vaporCoeffs[1]*pow(theta, 1.5) + 
                              vaporCoeffs[2]*pow(theta, 3.0) + 
                              vaporCoeffs[3]*pow(theta, 3.5) + 
                              vaporCoeffs[4]*pow(theta, 4.0) + 
                              vaporCoeffs[5]*pow(theta, 7.5)));
}

void Iawps95::saturatedCurve() const
{
    int n = 370;
    double Tstart = 275.0;
    double Tend = 645.0;
    double dT = (Tend - Tstart)/n;

    std::vector<double> temperature = std::vector<double>(n);
    std::vector<double> densityS = std::vector<double>(n);
    std::vector<double> energy = std::vector<double>(n);
    std::vector<double> pressure = std::vector<double>(n);
    std::vector<double> soundSpeed = std::vector<double>(n);


    for (int i = 0; i < n; i++)
    {
        temperature[i] = Tstart + i*dT;
        densityS[i] = saturatedVaporDensity(temperature[i]);

        double liquidDensity = saturatedLiquidDensity(temperature[i]);

        energy[i] = e(densityS[i], temperature[i]);
        soundSpeed[i] = std::sqrt(w2(densityS[i], temperature[i]));

        pressure[i] = vaporPressure(temperature[i]);
    }

    std::ofstream f;
	f.open("saturatedCurve.dat", std::ios::out);
	
	f << "#T, rho, e\n";

    for (int i = 0; i < temperature.size(); i++)
    {
        //f << temperature[i] << " " << densityS[i] << " " << 1.0/densityS[i] << " " << energy[i] << " " << pressure[i] << "\n";
        if(densityS[i] <= 125)
            f << densityS[i] << " " << energy[i] << " " << temperature[i] << " " << pressure[i] << " " << soundSpeed[i] << "\n";
    }

    f << std::endl;
    f.close();    
}

// dimensionless Helmholtz free energy functions

double Iawps95::phi0(double delta, double tau) const
{
    double out = log(delta) + coeffs.n0[0] + coeffs.n0[1]*tau + coeffs.n0[2]*log(tau);
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*log(1.0 - exp(-coeffs.gamma0[i]*tau));
    }
    return out;
}

double Iawps95::phi0d(double delta, double tau) const
{
    return 1.0/delta;
}

double Iawps95::phi0dd(double delta, double tau) const
{
    return -1.0/(delta*delta);
}

double Iawps95::phi0t(double delta, double tau) const
{
    double out = coeffs.n0[1] + coeffs.n0[2]/tau;
    for (int i = 3; i < 8; i++)
    {
        out += coeffs.n0[i]*coeffs.gamma0[i]*(1/(1.0 - exp(-coeffs.gamma0[i]*tau)) - 1.0);
    }    
    return out;
}

double Iawps95::phi0tt(double delta, double tau) const
{
    double out = -coeffs.n0[2]/(tau*tau);
    for (int i = 3; i < 8; i++)
    {
        double tmp = exp(-coeffs.gamma0[i]*tau);
        out -= coeffs.n0[i]*coeffs.gamma0[i]*coeffs.gamma0[i]*tmp*pow(1.0 - tmp, -2.0);
    }    
    return out;
}

double Iawps95::phi0dt(double delta, double tau) const
{
    return 0.0;
}

double Iawps95::phir(double delta, double tau) const
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

double Iawps95::phird(double delta, double tau) const
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

double Iawps95::phirdd(double delta, double tau) const
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

double Iawps95::phirt(double delta, double tau) const
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

double Iawps95::phirtt(double delta, double tau) const
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

double Iawps95::phirdt(double delta, double tau) const
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

double Iawps95::deltaFunc(double delta, double tau, int i) const
{
    return pow(thetaFunc(delta, tau, i), 2) + coeffs.B[i]*pow(pow(delta - 1.0, 2), coeffs.a[i]);
}

double Iawps95::thetaFunc(double delta, double tau, int i) const
{
    return 1.0 - tau + coeffs.A[i]*pow(pow(delta - 1.0, 2), 1.0/(2.0*coeffs.beta[i]));
}

double Iawps95::psiFunc(double delta, double tau, int i) const
{
    return exp(-coeffs.C[i]*pow(delta - 1.0, 2) - coeffs.D[i]*pow(tau - 1.0, 2));
}


double Iawps95::deltaFuncd(double delta, double tau, int i) const
{
    return (delta - 1.0)*(coeffs.A[i]*thetaFunc(delta, tau, i)*(2.0/coeffs.beta[i])*pow(pow(delta - 1.0, 2), 1/(2.0*coeffs.beta[i]) - 1.0) + 2.0*coeffs.B[i]*coeffs.a[i]*pow(pow(delta - 1.0, 2), coeffs.a[i] - 1.0));//pow(delta - 1.0, coeffs.a[i] - 1.0)
}

double Iawps95::deltaFuncdd(double delta, double tau, int i) const
{
    return (1.0/(delta - 1.0))*deltaFuncd(delta, tau, i) + pow(delta - 1.0, 2)*(4.0*coeffs.B[i]*coeffs.a[i]*(coeffs.a[i] - 1.0)*pow(pow(delta - 1.0, 2), coeffs.a[i] - 2.0) + 2.0*pow(coeffs.A[i], 2)*pow(1/coeffs.beta[i], 2)*pow(pow(pow(delta - 1.0, 2), 1/(2*coeffs.beta[i]) - 1.0), 2) + coeffs.A[i]*thetaFunc(delta, tau, i)*(4.0/coeffs.beta[i])*(1/(2*coeffs.beta[i]) - 1.0)*pow(pow(delta - 1.0, 2), 1/(2*coeffs.beta[i] - 2.0)));
}


double Iawps95::deltaFuncbid(double delta, double tau, int i) const
{
    return coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*deltaFuncd(delta, tau, i);
}

double Iawps95::deltaFuncbidd(double delta, double tau, int i) const
{
    return coeffs.b[i]*(pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*deltaFuncdd(delta, tau, i) + (coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0)*pow(deltaFuncd(delta, tau, i), 2));
}

double Iawps95::deltaFuncbit(double delta, double tau, int i) const
{
    return -2.0*thetaFunc(delta, tau, i)*coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0);
}

double Iawps95::deltaFuncbitt(double delta, double tau, int i) const
{
    return 2.0*coeffs.b[i]*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0) + 4.0*pow(thetaFunc(delta, tau, i), 2.0)*coeffs.b[i]*(coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0);
}

double Iawps95::deltaFuncbidt(double delta, double tau, int i) const
{
    return -coeffs.A[i]*coeffs.b[i]*(2.0/coeffs.beta[i])*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 1.0)*(delta - 1.0)*pow(pow(delta - 1, 2), 1.0/(2.0*coeffs.beta[i]) - 1.0) - 2.0*thetaFunc(delta, tau, i)*coeffs.b[i]*(coeffs.b[i] - 1.0)*pow(deltaFunc(delta, tau, i), coeffs.b[i] - 2.0)*deltaFuncd(delta, tau, i);
}


double Iawps95::psiFuncd(double delta, double tau, int i) const
{
    return -2.0*coeffs.C[i]*(delta - 1.0)*psiFunc(delta, tau, i);
}

double Iawps95::psiFuncdd(double delta, double tau, int i) const
{
    return (2.0*coeffs.C[i]*pow(delta - 1.0, 2) - 1.0)*2.0*coeffs.C[i]*psiFunc(delta, tau, i);
}

double Iawps95::psiFunct(double delta, double tau, int i) const
{
    return -2.0*coeffs.D[i]*(tau - 1.0)*psiFunc(delta, tau, i);
}

double Iawps95::psiFunctt(double delta, double tau, int i) const
{
    return (2.0*coeffs.D[i]*(tau - 1.0)*(tau - 1.0) - 1.0)*2.0*coeffs.D[i]*psiFunc(delta, tau, i);
}

double Iawps95::psiFuncdt(double delta, double tau, int i) const
{
    return 4.0*coeffs.C[i]*coeffs.D[i]*(delta - 1.0)*(tau - 1.0)*psiFunc(delta, tau, i);
}
