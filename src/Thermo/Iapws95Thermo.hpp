#ifndef IAPWS95THERMO
#define IAPWS95THERMO

#include "Thermo.hpp"
#include "StateEquations/Iapws95.hpp"

class Iapws95Thermo : public Thermo, public Iapws95
{
    public:

        Iapws95Thermo() : Thermo(), Iapws95() {}

        virtual ~Iapws95Thermo() {}

        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        //Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const;
        //Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const;

        //Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

        double calcM2is(double p0, double T0, double p2) const;
        
};

#endif // IAPWS95THERMO