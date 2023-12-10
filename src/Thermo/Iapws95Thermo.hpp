#ifndef IAPWS95THERMO
#define IAPWS95THERMO

#include "Thermo.hpp"
#include "StateEquations/Iapws95.hpp"

class Iapws95Thermo : public Thermo, Iapws95
{
    public:

        Iapws95Thermo() : Thermo(), Iapws95() {}

        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        virtual Vars<3> updateThermo(const Primitive& data) const = 0;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;
        Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;        
};

#endif // IAPWS95THERMO