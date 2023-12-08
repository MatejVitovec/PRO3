#include "FreeBoundary.hpp"

Compressible FreeBoundary::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    return w;
}