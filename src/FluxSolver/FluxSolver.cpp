#include <cmath>

#include "FluxSolver.hpp"

Field<Vars<5>> FluxSolver::claculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(wl.size());

    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], vector3toVars(faceList[i].normalVector))*faceList[i].area;
    }
    
    return out;
}