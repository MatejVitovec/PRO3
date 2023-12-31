#ifndef SPECIALGASTHERMO
#define SPECIALGASTHERMO

#include "Thermo.hpp"
#include "StateEquations/SpecialGas.hpp"

class SpecialGasThermo : public Thermo, SpecialGas
{
    public:

        SpecialGasThermo() : Thermo(), SpecialGas() {}

        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Vars<3> updateThermo(const Primitive& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;
        Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;        
};

#endif // SPECIALGASTHERMO