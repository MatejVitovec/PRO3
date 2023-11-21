#include <cmath>

#include "Limiter.hpp"

Field<Vars<5>> Limiter::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    return Field<Vars<5>>(mesh.getCellsSize(), Vars<5>());
    //return Field<Vars<5>>(mesh.getCellsSize(), Vars<5>({1.0, 1.0, 1.0, 1.0, 1.0}));
}