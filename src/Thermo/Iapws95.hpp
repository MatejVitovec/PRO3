#ifndef IAPWS95
#define IAPWS95

#include "Thermo.hpp"

class Iapws95 : public Thermo
{
    struct Coeffs
    {
        std::array<double, 8> n0;
        std::array<double, 8> gamma0;
        std::array<double, 56> c;
        std::array<double, 56> d;
        std::array<double, 56> t;
        std::array<double, 56> n;
        std::array<double, 56> alpha;
        std::array<double, 56> beta;
        std::array<double, 56> gamma;
        std::array<double, 56> epsilon;
        std::array<double, 56> a;
        std::array<double, 56> b;
        std::array<double, 56> B;
        std::array<double, 56> C;
        std::array<double, 56> D;
        std::array<double, 56> A;
    };

    public:

        Iapws95();

        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;
        Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

        double p(double rho, double T) const;
        double e(double rho, double T) const;
        double s(double rho, double T) const;
        double h(double rho, double T) const;
        double a2(double rho, double T) const;
        double a(double rho, double T) const;

        double tFromRhoE(double rho, double e, double guessT) const;
        double tFromRhoP(double rho, double p, double guessT) const;
        //double tFromRhoS(double rho, double s, double guessT) const;
        double rhoFromTP(double T, double p, double guessRho) const; //TODO

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;


    private:
        static constexpr double critT = 647.096;
        static constexpr double critRho = 322.0;
        static constexpr double specGasConst = 461.51805;

        static constexpr double numericalTolerance = 0.0001;

        Coeffs coeffs;

        std::vector<double> loadCoeffFile(std::string name, std::string dirName, int size) const;
        void loadCoeffs(std::string dirPath);

        double phi0(double delta, double tau) const;
        double phi0d(double delta, double tau) const;
        double phi0dd(double delta, double tau) const;
        double phi0t(double delta, double tau) const;
        double phi0tt(double delta, double tau) const;
        double phi0dt(double delta, double tau) const;

        double phir(double delta, double tau) const;
        double phird(double delta, double tau) const;
        double phirdd(double delta, double tau) const;
        double phirt(double delta, double tau) const;
        double phirtt(double delta, double tau) const;
        double phirdt(double delta, double tau) const;

        double deltaFunc(double delta, double tau, int i) const;
        double thetaFunc(double delta, double tau, int i) const;
        double psiFunc(double delta, double tau, int i) const;

        double deltaFuncd(double delta, double tau, int i) const;
        double deltaFuncdd(double delta, double tau, int i) const;

        double deltaFuncbid(double delta, double tau, int i) const;
        double deltaFuncbidd(double delta, double tau, int i) const;
        double deltaFuncbit(double delta, double tau, int i) const;
        double deltaFuncbitt(double delta, double tau, int i) const;
        double deltaFuncbidt(double delta, double tau, int i) const;

        double psiFuncd(double delta, double tau, int i) const;
        double psiFuncdd(double delta, double tau, int i) const;
        double psiFunct(double delta, double tau, int i) const;
        double psiFunctt(double delta, double tau, int i) const;
        double psiFuncdt(double delta, double tau, int i) const;
};

#endif // IAPWS95