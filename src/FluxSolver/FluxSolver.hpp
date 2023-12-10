#ifndef FLUXSOLVER_HPP
#define FLUXSOLVER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../Primitive.hpp"

class FluxSolver
{
    public:

        FluxSolver() {}

        virtual ~FluxSolver() {}

        Field<Vars<5>> calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const std::vector<Face>& faceList) const;
        Field<Vars<5>> calculateFluxes(const Field<Primitive>& ul, const Field<Primitive>& ur, const std::vector<Face>& faceList) const;
        virtual Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const = 0;
        virtual Vars<5> claculateFlux(const Primitive& ul, const Primitive& ur, const Vars<3>& normalVector) const = 0;

    protected:

};

#endif // FLUXSOLVER_HPP