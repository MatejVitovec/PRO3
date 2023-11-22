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

Field<Compressible> Thermo::updateInetrnalFieldFaces(Field<Compressible> wn, const Field<Compressible>& w, const Mesh& mesh) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            wn[i].setThermoVar(updateThermo(wn[i], w[i]));
        }
    }

    return wn;
}