#ifndef IAPWS95INTERPOLATIONTHERMO
#define IAPWS95INTERPOLATIONTHERMO

#include <memory>

#include "Thermo.hpp"
#include "StateEquations/Iapws95.hpp"
#include "Interpolation/Interpolation.hpp"

class Iapws95InterpolationThermo : public Thermo, Iapws95
{
    public:
        enum InterpolationType{BILINEAR, BIQUADRATIC};

        Iapws95InterpolationThermo(InterpolationType interpolationType);

        Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const;
        Vars<3> updateThermo(const Compressible& data) const;
        Vars<3> updateThermo(const Primitive& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;
        Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const;

    private:
        std::unique_ptr<Interpolation> pressureInterpolationFromRhoE;
        std::unique_ptr<Interpolation> soundSpeedInterpolationFromRhoE;
        std::unique_ptr<Interpolation> temperatureInterpolationFromRhoE;
        //std::unique_ptr<Interpolation> energyFromRhoPAux;
};

#endif // IAPWS95INTERPOLATIONTHERMO