#include "PressureTemperatureInlet.hpp"

Compressible PressureTemperatureInlet::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w);
    //return thermoModel->isentropicInletPressureTemperature(totalPressure, totalTemperature, velocityDirection, w); //TODO
}

void PressureTemperatureInlet::init(const Thermo * const thermoModel)
{
    std::array<double, 3> inletState = thermoModel->initPressureTemperatureInlet(totalPressure, totalTemperature);

    totalDensity = inletState[0];
    totalEntropy = inletState[1];
    totalEnthalpy = inletState[2];
}