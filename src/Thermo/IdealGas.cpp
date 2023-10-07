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

double IdealGas::pressure(const Compressible& data) const
{
    return (gamma - 1.0)*data.density()*(data.totalEnergy() - 0.5*data.absVelocity2());
}

double IdealGas::internalEnergy(const Compressible& data) const
{
    return pressure(data)/((gamma - 1.0)*data.density());
}

double IdealGas::soundSpeed(const Compressible& data) const
{
    return std::sqrt((gamma*pressure(data))/data.density());
}

double IdealGas::temperature(const Compressible& data) const
{
    return pressure(data)/(R*data[Compressible::RHO]);
}

Vars<3> IdealGas::calculateThermo(const Compressible& data) const
{
    return Vars<3>({temperature(data), pressure(data), soundSpeed(data)});
}

Compressible IdealGas::primitiveToConservative(const Vars<5>& primitive) const
{
    double velocity2 = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

    return Compressible({primitive[0],
                         primitive[0]*primitive[1],
                         primitive[0]*primitive[2],
                         primitive[0]*primitive[3],
                         0.5*primitive[0]*velocity2 + (primitive[4])/(gamma - 1.0)});
}