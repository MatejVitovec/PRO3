#ifndef SPECIALGASTHERMO
#define SPECIALGASTHERMO

#include "Thermo.hpp"
#include "StateEquations/SpecialGas.hpp"

class SpecialGasThermo : public Thermo, SpecialGas
{
    public:

        SpecialGasThermo() : Thermo(), SpecialGas() {}

        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        //Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const;
        //Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const;

        //Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;
};

#endif // SPECIALGASTHERMO