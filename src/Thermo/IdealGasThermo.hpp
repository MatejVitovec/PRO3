#ifndef IDEALGASTHERMO_HPP
#define IDEALGASTHERMO_HPP

#include "Thermo.hpp"

class IdealGasThermo : public Thermo
{
    public:

        IdealGasThermo() : Thermo(), gamma(1.4), R(287.05) {}
        IdealGasThermo(double gamma_, double R_) : Thermo(), gamma(gamma_), R(R_) {}

        //overwritten virtual
        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;
        Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

    private:
        double gamma;
        double R;

        double pressure(const Compressible& data) const;
        double internalEnergy(const Compressible& data) const;
        double soundSpeed(const Compressible& data) const;
        double temperature(const Compressible& data) const;

};

#endif // IDEALGASTHERMO_HPP