#include "FreeBoundary.hpp"

Compressible FreeBoundary::calculateState(const Compressible& wl, const Compressible& wr, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = wl;

    out.setThermoVar(thermoModel->updateThermo(out, wr));

    return out;
}