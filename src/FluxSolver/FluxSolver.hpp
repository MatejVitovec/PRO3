#ifndef FLUXSOLVER_HPP
#define FLUXSOLVER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

class FluxSolver
{
    public:

        FluxSolver() {}

        virtual ~FluxSolver() {}

        Field<Vars<5>> claculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const std::vector<Face>& faceList) const;
        virtual Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const = 0;

    protected:

};

#endif // FLUXSOLVER_HPP