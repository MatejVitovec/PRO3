#include <cmath>

#include "Thermo.hpp"

Field<Compressible> Thermo::updateField(const Field<Compressible>& w) const
{
    Field<Compressible> out(w.size());

    for (int i = 0; i < w.size(); i++)
    {
        w[i].setThermoVar(calculateThermo(w[i]));
    }

    return out;
}