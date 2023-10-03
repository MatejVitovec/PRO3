#ifndef IDEALGAS_HPP
#define IDEALGAS_HPP

#include "Thermo.hpp"

class IdealGas : public Thermo
{
    public:

        IdealGas() : Thermo(), gamma(1.4), R(287.05) {}
        IdealGas(double gamma_, double R_, double cp_) : Thermo(), gamma(gamma_), R(R_) {}

        void setGamma(double gamma_);
        void setR(double R_);
        double getGamma() const;
        double getR() const;

        //overwritten virtual
        Vars<3> calculateThermo(const Compressible& data) const;

    private:
        double gamma;
        double R;

        double pressure(const Compressible& data) const;
        double internalEnergy(const Compressible& data) const;
        double soundSpeed(const Compressible& data) const;

        double density(double totalDensity, double machNumeber2) const;
        double soundSpeed(double pressure, double density) const;
        double machNumber2(double totalPressure, double pressure) const;
        double cp() const;

};

#endif // IDEALGAS_HPP