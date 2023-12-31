#ifndef HLL_HPP
#define HLL_HPP

#include "FluxSolver.hpp"

class Hll : public FluxSolver
{
    public:

        Hll() {}

        virtual ~Hll() {}

        Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const;
        Vars<5> claculateFlux(const Primitive& wl, const Primitive& wr, const Vars<3>& normalVector) const {return Vars<5>(); };

    private:
        Vars<3> waveSpeedsEstimate(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const;

};

#endif // HLL_HPP