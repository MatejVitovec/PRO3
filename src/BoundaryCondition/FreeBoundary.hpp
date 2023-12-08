#ifndef FREEBOUNDARY_HPP
#define FREEBOUNDARY_HPP

#include "BoundaryCondition.hpp"

class FreeBoundary : public BoundaryCondition
{
    public:

        FreeBoundary(Boundary meshBoundary) : BoundaryCondition(meshBoundary, FREEBOUNDARY) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;

};

#endif // FREEBOUNDARY