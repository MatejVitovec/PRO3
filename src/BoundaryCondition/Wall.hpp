#ifndef WALL_HPP
#define WALL_HPP

#include "BoundaryCondition.hpp"

class Wall : public BoundaryCondition
{
    public:

        Wall(Boundary meshBoundary) : BoundaryCondition(meshBoundary, WALL) {}

        Compressible calculateState(const Compressible& wl, const Compressible& wr, const Face& f, const Thermo * const thermoModel) const;
    
    private:

};

#endif // WALL