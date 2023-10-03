#include <cmath>

#include "IdealGas.hpp"

void IdealGas::setGamma(double gamma_)
{
    gamma = gamma_;
}

void IdealGas::setR(double R_)
{
    R = R_;
}

double IdealGas::getGamma() const
{
    return gamma;
}

double IdealGas::getR() const
{
    return R;
}

double IdealGas::cp() const
{
    return (gamma*R)/(gamma - 1.0);
}

double IdealGas::pressure(const Compressible& data) const
{
    return (gamma - 1.0)*data.density()*(data.totalEnergy() - 0.5*data.absVelocity2());
}

double IdealGas::internalEnergy(const Compressible& data) const
{
    return data.pressure()/((gamma - 1.0)*data.density());
}

double IdealGas::soundSpeed(const Compressible& data) const
{
    return std::sqrt((gamma*data.pressure())/data.density());
}

double IdealGas::density(double totalDensity, double machNumeber2) const
{
    return totalDensity*std::pow(1.0 + ((gamma - 1.0)/2.0)*machNumeber2, 1.0/(1.0 - gamma));
}

double IdealGas::soundSpeed(double pressure, double density) const
{
    return std::sqrt((gamma*pressure)/density);
}

double IdealGas::machNumber2(double totalPressure, double pressure) const
{
    return (2.0/(gamma - 1.0))*(std::pow((pressure/totalPressure), ((1.0 - gamma)/gamma)) - 1.0);
}

Vars<3> IdealGas::calculateThermo(const Compressible& data) const
{
    return Vars<3>({pressure(data), pressure(data), soundSpeed(data)});
}