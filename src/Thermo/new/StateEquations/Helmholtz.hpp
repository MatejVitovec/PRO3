#ifndef HELMHOLTZ
#define HELMHOLTZ

#include <array>
#include <vector>

#include "../../NonLinearSolver/NewtonMethod.hpp"

class Helmholtz
{
    public:

        Helmholtz();

        virtual ~Helmholtz() {}

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

        NewtonMethod nonLinearSolver;

        virtual double phi0(double delta, double tau) const = 0;
        virtual double phi0d(double delta, double tau) const = 0;
        virtual double phi0dd(double delta, double tau) const = 0;
        virtual double phi0t(double delta, double tau) const = 0;
        virtual double phi0tt(double delta, double tau) const = 0;
        virtual double phi0dt(double delta, double tau) const = 0;

        virtual double phir(double delta, double tau) const = 0;
        virtual double phird(double delta, double tau) const = 0;
        virtual double phirdd(double delta, double tau) const = 0;
        virtual double phirt(double delta, double tau) const = 0;
        virtual double phirtt(double delta, double tau) const = 0;
        virtual double phirdt(double delta, double tau) const = 0;

        std::function<double(double, double)> pFunc = [this](double delta, double tau) { return critT*critRho*specGasConst*(delta/tau)*(1.0 + delta*phird(delta, tau)); };
        std::function<double(double, double)> eFunc = [this](double delta, double tau) { return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau)); };
        std::function<double(double, double)> sFunc = [this](double delta, double tau) { return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau)); };

        std::function<double(double, double)> pDDeltaFunc = [this](double delta, double tau) { return critRho*critT*(specGasConst/tau)*(1.0 + 2.0*delta*phird(delta, tau) + delta*delta*phirdd(delta, tau)); };
        std::function<double(double, double)> sDDeltaFunc = [this](double delta, double tau) { return specGasConst*(tau*(phi0dt(delta, tau) + phirdt(delta, tau)) - phi0d(delta, tau) - phird(delta, tau)); };

        std::function<double(double, double)> pDTauFunc = [this](double delta, double tau) { return specGasConst*critRho*critT*(delta/tau)*((-1.0/tau) - (delta/tau)*phird(delta, tau) + delta*phirdt(delta, tau)); };
        std::function<double(double, double)> eDTauFunc = [this](double delta, double tau) { return specGasConst*critT*(phi0tt(delta, tau) + phirtt(delta, tau)); };
        std::function<double(double, double)> sDTauFunc = [this](double delta, double tau) { return specGasConst*tau*(phi0tt(delta, tau) + phirtt(delta, tau)); };

};

#endif // HELMHOLTZ