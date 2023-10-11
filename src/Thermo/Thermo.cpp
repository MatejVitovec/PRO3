#include <cmath>

#include "Thermo.hpp"

Field<Compressible> Thermo::updateField(Field<Compressible> wn, const Field<Compressible>& w) const
{
    for (int i = 0; i < wn.size(); i++)
    {
        wn[i].setThermoVar(updateThermo(wn[i], w[i]));
    }

    return wn;
}