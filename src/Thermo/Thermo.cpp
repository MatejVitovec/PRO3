#include <cmath>

#include "Thermo.hpp"

Field<Compressible> Thermo::updateField(Field<Compressible> w) const
{
    for (int i = 0; i < w.size(); i++)
    {
        w[i].setThermoVar(updateThermo(w[i]));
    }

    return w;
}