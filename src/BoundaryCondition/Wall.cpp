#include "Wall.hpp"

Compressible Wall::calculateState(const Compressible& wl, const Compressible& wrOld, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = wl;

    Vars<3> normalVector = vector3toVars(f.normalVector);
    Vars<3> ghostVelocity = wl.velocity() - 2*wl.normalVelocity(normalVector)*normalVector;
    double density = wl.density();

    out[Compressible::RHO_U] = density*ghostVelocity[0];
    out[Compressible::RHO_V] = density*ghostVelocity[1];
    out[Compressible::RHO_W] = density*ghostVelocity[2];

    out.setThermoVar(thermoModel->updateThermo(out, wrOld));

    return out;
}

void Wall::correct(const Field<Compressible>& w, Field<Compressible> wl, Field<Compressible>& wr, const Field<Compressible>& wrOld, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (auto & faceIndex : boundary.facesIndex)
    {
        Vars<5> wlDiff = dot(grad[ownerIndexList[faceIndex]], vector3toVars(faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center));

        wl[faceIndex] = w[ownerIndexList[faceIndex]] + phi[ownerIndexList[faceIndex]]*wlDiff;

        wr[faceIndex] = calculateState(wl[faceIndex], wrOld[faceIndex], faces[faceIndex], thermoModel);
    }
}