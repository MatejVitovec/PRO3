#include "PressureTemperatureInlet.hpp"


void PressureTemperatureInlet::setTotalPressure(double totalPressure_)
{
    totalPressure = totalPressure_;
}

void PressureTemperatureInlet::setTotalTemperature(double totalTemperature_)
{
    totalTemperature = totalTemperature_;
}

void PressureTemperatureInlet::setVelocityDirection(Vars<3> velocityDirection_)
{
    velocityDirection = velocityDirection_;
}

double PressureTemperatureInlet::getTotalPressure() const
{
    return totalPressure;
}

double PressureTemperatureInlet::getTotalTemperature() const
{
    return totalTemperature;
}

Vars<3> PressureTemperatureInlet::getVelocityDirection() const
{
    return velocityDirection;
}

Compressible PressureTemperatureInlet::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, velocityDirection, wl);
}