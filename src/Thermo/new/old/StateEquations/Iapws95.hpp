#ifndef IAPWS95
#define IAPWS95

#include "Iapws95Base.hpp"

class Iapws95 : public Iapws95Base
{
    struct Coeffs
    {
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

    private:
    
        Coeffs coeffs;

        std::vector<double> loadCoeffFile(std::string name, std::string dirName, int size) const;
        void loadCoeffs(std::string dirPath);

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