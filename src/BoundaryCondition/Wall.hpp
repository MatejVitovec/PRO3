#ifndef WALL_HPP
#define WALL_HPP

#include "BoundaryCondition.hpp"

class Wall : public BoundaryCondition
{
    public:

        Wall(Boundary meshBoundary) : BoundaryCondition(meshBoundary, WALL) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;

        void correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;

    private:

};

#endif // WALL