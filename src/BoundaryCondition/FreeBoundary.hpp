#ifndef FREEBOUNDARY_HPP
#define FREEBOUNDARY_HPP

#include "BoundaryCondition.hpp"

class FreeBoundary : public BoundaryCondition
{
    public:

        FreeBoundary() {}
        FreeBoundary(Boundary meshBoundary) : BoundaryCondition(meshBoundary) {}

        Compressible calculateState(const Compressible& wl, const Face& f) const;

};

#endif // FREEBOUNDARY