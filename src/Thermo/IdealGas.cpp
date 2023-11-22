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
    return std::sqrt((gamma*std::max(0.0, pressure(data)))/data.density());
}

double IdealGas::temperature(const Compressible& data) const
{
    return pressure(data)/(R*data[Compressible::RHO]);
}

Vars<3> IdealGas::updateThermo(const Compressible& data, const Compressible& dataOld) const
{
    return Vars<3>({temperature(data), pressure(data), soundSpeed(data)});
}

Vars<3> IdealGas::updateThermo(const Compressible& data) const
{
    return Vars<3>({temperature(data), pressure(data), soundSpeed(data)});
}

Compressible IdealGas::primitiveToConservative(const Vars<5>& primitive) const
{
    double velocity2 = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

    Compressible out = Compressible({primitive[0],
                                     primitive[0]*primitive[1],
                                     primitive[0]*primitive[2],
                                     primitive[0]*primitive[3],
                                     0.5*primitive[0]*velocity2 + (primitive[4])/(gamma - 1.0)});

    out.setThermoVar(Vars<3>({temperature(out), pressure(out), soundSpeed(out)}));

    return out;
}

Compressible IdealGas::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const
{
    double p = std::min(pressure(stateIn), pTot);
    double M2 = (2.0/(gamma - 1.0))*(std::pow((pTot/p), ((gamma - 1.0)/gamma)) - 1.0);
    double T = TTot/(1.0 + ((gamma - 1.0)/2)*M2);
    double rho = p/(R*T);
    double a = std::sqrt((gamma*p)/rho);
    double absU= std::sqrt(M2)*a;

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         0.5*rho*absU*absU + p/(gamma - 1.0)},
                         {T, p, a});
}

Compressible IdealGas::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const
{
    //TODO
    return Compressible();
}