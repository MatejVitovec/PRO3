#include "Wall.hpp"

Compressible Wall::calculateState(const Compressible& wl, const Face& f) const
{
    Vars<3> normalVector = vector3toVars(f.normalVector);
    Vars<3> ghostVelocity = wl.velocity() - 2*wl.normalVelocity(normalVector)*normalVector;
    double density = wl.density();

    return Compressible({density, density*ghostVelocity[0], density*ghostVelocity[1], density*ghostVelocity[2], density*wl.totalEnergy()});
}