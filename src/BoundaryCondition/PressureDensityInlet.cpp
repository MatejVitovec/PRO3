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

Compressible PressureDensityInlet::calculateState(const Compressible& wl, const Face& f) const
{
    std::shared_ptr<EquationOfState> eqs = Compressible::getEquationOfState();

    double pressure = std::min(wl.pressure(), totalPressure); //interpolate from domain

    double machNumber2 = eqs->machNumber2(totalPressure, pressure);
    double density = eqs->density(totalDensity, machNumber2);
    double absVelocity = std::sqrt(machNumber2)*eqs->soundSpeed(pressure, density);

    return eqs->primitiveToConservative(Vars<5>({density,
                                                 absVelocity*velocityDirection[0],
                                                 absVelocity*velocityDirection[1],
                                                 absVelocity*velocityDirection[2],
                                                 pressure}));
}