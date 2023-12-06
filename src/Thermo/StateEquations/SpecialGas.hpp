#ifndef SPECIALGAS
#define SPECIALGAS

#include "Helmholtz.hpp"

class SpecialGas : public Helmholtz
{
    public:

        SpecialGas();

        virtual ~SpecialGas() {}

    private:

    struct Coeffs
    {
        std::array<double, 8> n0;
        std::array<double, 8> gamma0;
        std::array<double, 7> d;
        std::array<double, 7> t;
        std::array<double, 7> n;
    };
    
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

#endif // SPECIALGAS