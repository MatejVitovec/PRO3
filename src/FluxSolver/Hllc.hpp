#ifndef HLLC_HPP
#define HLLC_HPP

#include "FluxSolver.hpp"

class Hllc : public FluxSolver
{
    public:

        Hllc() {}

        virtual ~Hllc() {}

        Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const;

};

#endif // HLLC_HPP