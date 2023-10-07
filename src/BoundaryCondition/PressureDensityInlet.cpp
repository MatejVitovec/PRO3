#include "PressureDensityInlet.hpp"


void PressureDensityInlet::setTotalPressure(double totalPressure_)
{
    totalPressure = totalPressure_;
}

void PressureDensityInlet::setTotalDensity(double totaltemperature_)
{
    totalDensity = totaltemperature_;
}

void PressureDensityInlet::setVelocityDirection(Vars<3> velocityDirection_)
{
    velocityDirection = velocityDirection_;
}

double PressureDensityInlet::getTotalPressure() const
{
    return totalPressure;
}

double PressureDensityInlet::getTotalDensity() const
{
    return totalDensity;
}

Vars<3> PressureDensityInlet::getVelocityDirection() const
{
    return velocityDirection;
}

Compressible PressureDensityInlet::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    double pressure = std::min(wl.pressure(), totalPressure); //interpolate from domain

    //TODO
    double gamma = 1.4;
    double machNumber2 = (2.0/(gamma - 1.0))*(std::pow((pressure/totalPressure), ((1.0 - gamma)/gamma)) - 1.0);
    double density = totalDensity*std::pow(1.0 + ((gamma - 1.0)/2.0)*machNumber2, 1.0/(1.0 - gamma));
    double absVelocity = std::sqrt(machNumber2)*std::sqrt((gamma*pressure)/density);

    return thermoModel->primitiveToConservative(Vars<5>({density,
                                                         absVelocity*velocityDirection[0],
                                                         absVelocity*velocityDirection[1],
                                                         absVelocity*velocityDirection[2],
                                                         pressure}));
}