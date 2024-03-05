#include <cmath>
#include <omp.h>

#include "Thermo.hpp"

Field<Compressible> Thermo::updateField(Field<Compressible> w) const
{
    #pragma omp parallel for
    for (int i = 0; i < w.size(); i++)
    {
        w[i].setThermoVar(updateThermo(w[i]));
    }

    return w;
}


/*Field<Compressible> Thermo::updateField(Field<Compressible> wn, const Field<Compressible>& w) const
{
    #pragma omp parallel for
    for (int i = 0; i < wn.size(); i++)
    {
        wn[i].setThermoVar(updateThermo(wn[i], w[i]));
    }

    return wn;
}*/

Field<Compressible> Thermo::updateInetrnalFieldFaces(Field<Compressible> w, const Mesh& mesh) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    #pragma omp parallel for
    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            w[i].setThermoVar(updateThermo(w[i]));
        }
    }

    return w;
}

/*Field<Compressible> Thermo::updateInetrnalFieldFaces(Field<Compressible> wn, const Field<Compressible>& w, const Mesh& mesh) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    #pragma omp parallel for
    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            wn[i].setThermoVar(updateThermo(wn[i], w[i]));
        }
    }

    return wn;
}*/


Field<Primitive> Thermo::updateField(Field<Primitive> u) const
{
    #pragma omp parallel for
    for (int i = 0; i < u.size(); i++)
    {
        u[i].setThermo(updateThermo(u[i]));
    }

    return u;
}

Field<Primitive> Thermo::updateInetrnalFieldFaces(Field<Primitive> u, const Mesh& mesh) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    #pragma omp parallel for
    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            u[i].setThermo(updateThermo(u[i]));
        }
    }

    return u;
}