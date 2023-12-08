#include "PressureTemperatureInlet.hpp"

Compressible PressureTemperatureInlet::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInletPressureTemperature(totalPressure, totalTemperature, velocityDirection, w, w); //TODO
}