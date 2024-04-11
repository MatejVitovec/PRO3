#include "IsentropicInlet.hpp"

Compressible IsentropicInlet::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w);
}

