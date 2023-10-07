#include "FreeBoundary.hpp"

Compressible FreeBoundary::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    return wl;
}