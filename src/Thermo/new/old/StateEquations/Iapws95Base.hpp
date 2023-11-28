#ifndef IAPWS95BASE
#define IAPWS95BASE

#include <array>
#include <vector>

#include "../../NonLinearSolver/NewtonMethod.hpp"

class Iapws95Base
{
    public:

        Iapws95Base();

        struct CoeffsIdealGasPart
        {
            std::array<double, 8> n0;
            std::array<double, 8> gamma0;
        };

        double p(double rho, double T) const;
        double e(double rho, double T) const;
        double s(double rho, double T) const;
        double h(double rho, double T) const;
        double a2(double rho, double T) const;
        double a(double rho, double T) const;

        double tFromRhoE(double rho, double e, double guessT) const;
        double tFromRhoP(double rho, double p, double guessT) const;
        double rhoFromTP(double T, double p, double guessRho) const;

        std::pair<double, double> rRhoFromSP(double s, double p, double guessRho, double guessT) const;

    protected:
        static constexpr double critT = 647.096;
        static constexpr double critRho = 322.0;
        static constexpr double specGasConst = 461.51805;

        CoeffsIdealGasPart coeffs0;

        NewtonMethod nonLinearSolver;

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

        std::function<double(double, double)> pFunc = [this](double delta, double tau) { return delta*critRho*specGasConst*(critT/tau)*(1.0 + delta*phird(delta, tau)); };
        std::function<double(double, double)> eFunc = [this](double delta, double tau) { return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau)); };
        std::function<double(double, double)> sFunc = [this](double delta, double tau) { return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau)); };

        std::function<double(double, double)> pDTauFunc = [this](double delta, double tau) { return specGasConst*critRho*critT*(delta/tau)*(1.0/tau + (delta/tau)*specGasConst*critRho*critT*(delta/tau) - delta*phirdt(delta, tau)); };        
        std::function<double(double, double)> eDTauFunc = [this](double delta, double tau) { return specGasConst*critT*(phi0tt(delta, tau) + phirtt(delta, tau)); };
        std::function<double(double, double)> sDTauFunc = [this](double delta, double tau) { return specGasConst*tau*(phi0tt(delta, tau) + phirtt(delta, tau)); };
        
        std::function<double(double, double)> pDDeltaFunc = [this](double delta, double tau) { return critRho*critT*(specGasConst/tau)*(1.0 + 2.0*delta*phird(delta, tau) + delta*delta*phirdt(delta, tau)); };
        std::function<double(double, double)> sDDeltaFunc = [this](double delta, double tau) { return specGasConst*(tau*(phi0dt(delta, tau) + phirdt(delta, tau)) - phi0d(delta, tau) - phird(delta, tau)); };

};

#endif // IAPWS95BASE