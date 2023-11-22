#include "PressureDensityInlet.hpp"


void PressureDensityInlet::setTotalPressure(double totalPressure_)
{
    totalPressure = totalPressure_;
}

void PressureDensityInlet::setTotalDensity(double totalDensity)
{
    totalDensity = totalDensity;
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

Compressible PressureDensityInlet::calculateState(const Compressible& wl, const Compressible& wr, const Face& f, const Thermo * const thermoModel) const
{
    //TODO

    double pressure = std::min(wl.pressure(), totalPressure); //interpolate from domain
    double gamma = 1.4;
    double machNumber2 = (2.0/(gamma - 1.0))*(std::pow((pressure/totalPressure), ((1.0 - gamma)/gamma)) - 1.0);
    //double machNumber2 = (2.0/(gamma - 1.0))*(std::pow((totalPressure/pressure), ((gamma - 1.0)/gamma)) - 1.0);
    double density = totalDensity*std::pow(1.0 + ((gamma - 1.0)/2.0)*machNumber2, 1.0/(1.0 - gamma));
    double absVelocity = std::sqrt(machNumber2)*std::sqrt((gamma*pressure)/density);

    return thermoModel->primitiveToConservative(Vars<5>({density,
                                                         absVelocity*velocityDirection[0],
                                                         absVelocity*velocityDirection[1],
                                                         absVelocity*velocityDirection[2],
                                                         pressure}));
}