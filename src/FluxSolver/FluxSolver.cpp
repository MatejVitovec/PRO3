#include <cmath>

#include "FluxSolver.hpp"

Field<Vars<5>> FluxSolver::calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(wl.size());

    #pragma omp parallel for
    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], vector3toVars(faceList[i].normalVector))*faceList[i].area;
    }
    
    return out;
}

Field<Vars<5>> FluxSolver::calculateFluxes(const Field<Primitive>& ul, const Field<Primitive>& ur, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(ul.size());

    #pragma omp parallel for
    for (int i = 0; i < ul.size(); i++)
    {
        out[i] = claculateFlux(ul[i], ur[i], vector3toVars(faceList[i].normalVector))*faceList[i].area;
    }
    
    return out;
}