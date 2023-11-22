#include "PressureOutlet.hpp"


void PressureOutlet::setPressure(double pressure_)
{
    pressure = pressure_;
}

double PressureOutlet::getPressure() const
{
    return pressure;
}

Compressible PressureOutlet::calculateState(const Compressible& wl, const Compressible& wr, const Face& f, const Thermo * const thermoModel) const
{
    Compressible aux = wl;
    aux.setThermoVar(thermoModel->updateThermo(aux, wr));

    if(wl.absVelocity()/aux.soundSpeed() >= 1.0)
    {
        return aux;
    }

    return thermoModel->primitiveToConservative(Vars<5>({wl.density(),
                                                         wl.velocityU(),
                                                         wl.velocityV(),
                                                         wl.velocityW(),
                                                         pressure}));
}