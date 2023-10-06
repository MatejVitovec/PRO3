#ifndef IAWPS95
#define IAWPS95

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

class Iawps95
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

        Iawps95();

        double p(double rho, double T) const;
        double e(double rho, double T) const;
        double s(double rho, double T) const;
        double h(double rho, double T) const; //TODO
        double w2(double rho, double T) const; //TODO

        double temperatureFromRhoE(double rho, double e, double guessT) const;
        double temperatureFromRhoP(double rho, double p, double guessT) const;

        void test(double rho, double T) const;



    private:
        static constexpr  double critT = 647.096;
        static constexpr  double critRho = 322.0;
        static constexpr  double specGasConst = 0.46151805;

        static constexpr  double numericalTolerance = 0.0001;

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

#endif // IAWPS95