#ifndef IDEALGASTHERMO_HPP
#define IDEALGASTHERMO_HPP

#include "Thermo.hpp"

class IdealGasThermo : public Thermo
{
    public:

        IdealGasThermo() : Thermo(), gamma(1.3326), R(461.51805), energyShift(1995917.326) {}
        IdealGasThermo(double gamma_, double R_) : Thermo(), gamma(gamma_), R(R_), energyShift(1995917.326) {}

        //overwritten virtual
        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        //Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const;
       // Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

    private:
        double gamma;
        double R;
        double energyShift;

        double pressure(const Compressible& data) const;
        double soundSpeed(const Compressible& data) const;
        double temperature(const Compressible& data) const;

};

#endif // IDEALGASTHERMO_HPP