#include "PressureOutlet.hpp"


Compressible PressureOutlet::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    if(w.absVelocity()/w.soundSpeed() >= 1.0)
    {
        return w;
    }

    return thermoModel->primitiveToConservative(Vars<5>({w.density(),
                                                         w.velocityU(),
                                                         w.velocityV(),
                                                         w.velocityW(),
                                                         pressure}));
}