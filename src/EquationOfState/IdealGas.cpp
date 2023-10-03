#include <cmath>

#include "IdealGas.hpp"

//musi byt pod IdealGas.cpp kvuli forward deklaraci - byl pouzit prototyp
#include "../Compressible.hpp"

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

Compressible IdealGas::primitiveToConservative(const Vars<5>& primitive) const
{
    double velocity2 = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

    return Compressible({primitive[0],
                         primitive[0]*primitive[1],
                         primitive[0]*primitive[2],
                         primitive[0]*primitive[3],
                         0.5*primitive[0]*velocity2 + (primitive[4])/(gamma - 1.0)});
}

Compressible IdealGas::nonLinearSubsonicInletBoundaryState(const Compressible& inDomainState, const Vars<3>& normalVector, double totalPressure, double totaltemperature, Vars<3> velocityDirection) const
{
    double gamma1 = gamma - 1.0;

    double inDomainSS = inDomainState.soundSpeed();
    Vars<3> inDomainVelocity = inDomainState.velocity();
    double riemannInv = dot(inDomainVelocity, normalVector) - (2*inDomainSS)/gamma1;
    
    double absVelo = inDomainState.absVelocity();
    
    //double cosPhi = inDomainState.velocityCosine(normalVector);
    double cosPhi = -((dot(inDomainVelocity, normalVector))/absVelo);
    if (absVelo == 0.0) cosPhi = 1.0;

    double stagnationSS2 = inDomainSS*inDomainSS + (gamma1/2)*inDomainState.absVelocity2();

    double soundSpeedBoundary = (((-riemannInv)*gamma1)/(gamma1*cosPhi*cosPhi + 2.0))*(1.0 + cosPhi*std::sqrt((((gamma1*cosPhi*cosPhi + 2.0)*stagnationSS2)/(gamma1*(riemannInv*riemannInv))) - gamma1/2.0));

    double boundaryTemperature = totaltemperature*((soundSpeedBoundary*soundSpeedBoundary)/stagnationSS2);

    double boundaryPressure = totalPressure*std::pow(boundaryTemperature*totaltemperature, gamma/gamma1);
    double boundarydensity = boundaryPressure/(R*boundaryTemperature);
    double boundaryAbsVelocity = std::sqrt(2*cp()*(totaltemperature - boundaryTemperature));

    return Compressible(primitiveToConservative(Vars<5>({boundarydensity,
                                                         boundaryAbsVelocity*velocityDirection[0],
                                                         boundaryAbsVelocity*velocityDirection[1],
                                                         boundaryAbsVelocity*velocityDirection[2],
                                                         boundaryPressure})));
}