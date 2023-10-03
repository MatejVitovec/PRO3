#include "PressureOutlet.hpp"


void PressureOutlet::setPressure(double pressure_)
{
    pressure = pressure_;
}

double PressureOutlet::getPressure() const
{
    return pressure;
}

Compressible PressureOutlet::calculateState(const Compressible& wl, const Face& f) const
{
    Vars<3> normalVector = vector3toVars(f.normalVector);
    double referenceSoundSpeed = wl.soundSpeed();

    if(wl.absVelocity()/referenceSoundSpeed >= 1.0)
    {
        return wl;
    }

    /*double referenceDensity = wl.density();
    double inDomainPressure = wl.pressure();

    double nonConstPressure = pressure; //pressure is const - non const cast copy

    double aux = (inDomainPressure - pressure)/(referenceDensity*referenceSoundSpeed);

    return Compressible::primitiveToConservative(Vars<5>({referenceDensity + (pressure - inDomainPressure)/(referenceSoundSpeed*referenceSoundSpeed),
                                                          wl.velocityU() + normalVector[0]*aux,
                                                          wl.velocityV() + normalVector[1]*aux,
                                                          wl.velocityW() + normalVector[2]*aux,
                                                          nonConstPressure}));*/

    return Compressible::primitiveToConservative(Vars<5>({wl.density(),
                                                          wl.velocityU(),
                                                          wl.velocityV(),
                                                          wl.velocityW(),
                                                          pressure}));
}