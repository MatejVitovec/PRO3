#include "Wall.hpp"

Compressible Wall::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = wl;

    Vars<3> normalVector = vector3toVars(f.normalVector);
    Vars<3> ghostVelocity = wl.velocity() - 2*wl.normalVelocity(normalVector)*normalVector;
    double density = wl.density();

    out[Compressible::RHO_U] = density*ghostVelocity[0];
    out[Compressible::RHO_V] = density*ghostVelocity[1];
    out[Compressible::RHO_W] = density*ghostVelocity[2];

    return out;
}