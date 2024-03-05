#ifndef IAPWS95SPECIALGAS
#define IAPWS95SPECIALGAS

#include "Thermo.hpp"

class Iapws95SpecialGas : public Thermo
{
    struct Coeffs
    {
        std::array<double, 8> n0;
        std::array<double, 8> gamma0;
        std::array<double, 7> d;
        std::array<double, 7> t;
        std::array<double, 7> n;
    };

    public:

        Iapws95SpecialGas();

        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTote) const;
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
};

#endif // IAPWS95SPECIALGAS