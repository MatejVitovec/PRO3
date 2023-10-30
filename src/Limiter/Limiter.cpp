#include <cmath>

#include "Limiter.hpp"

Field<Vars<5>> Limiter::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, Field<std::array<Vars<3>, 5>> grad, const Mesh& mesh) const
{
    return Field<Vars<5>>(mesh.getCellsSize(), Vars<5>({0.0, 0.0, 0.0, 0.0, 0.0}));
    //return Field<Vars<5>>(mesh.getCellsSize(), Vars<5>({1.0, 1.0, 1.0, 1.0, 1.0}));
}