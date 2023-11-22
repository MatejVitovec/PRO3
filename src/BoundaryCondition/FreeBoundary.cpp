#include "FreeBoundary.hpp"

Compressible FreeBoundary::calculateState(const Compressible& wl, const Compressible& wrOld, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = wl;

    out.setThermoVar(thermoModel->updateThermo(out, wrOld));

    return out;
}