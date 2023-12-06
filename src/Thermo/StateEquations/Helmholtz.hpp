#ifndef HELMHOLTZ
#define HELMHOLTZ

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

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

        std::pair<double, double> RhoTFromSP(double s, double p, double guessRho, double guessT) const;

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

        double pFunc(double delta, double tau) const;
        double eFunc(double delta, double tau) const;
        double sFunc(double delta, double tau) const;

        double pDDeltaFunc(double delta, double tau) const;
        double sDDeltaFunc(double delta, double tau) const;

        double pDTauFunc(double delta, double tau) const;
        double eDTauFunc(double delta, double tau) const;
        double sDTauFunc(double delta, double tau) const;

};

#endif // HELMHOLTZ