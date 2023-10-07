#include "PressureOutlet.hpp"


void PressureOutlet::setPressure(double pressure_)
{
    pressure = pressure_;
}

double PressureOutlet::getPressure() const
{
    return pressure;
}

Compressible PressureOutlet::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    double referenceSoundSpeed = wl.soundSpeed();

    if(wl.absVelocity()/referenceSoundSpeed >= 1.0)
    {
        return wl;
    }

    return thermoModel->primitiveToConservative(Vars<5>({wl.density(),
                                                         wl.velocityU(),
                                                         wl.velocityV(),
                                                         wl.velocityW(),
                                                         pressure}));
}