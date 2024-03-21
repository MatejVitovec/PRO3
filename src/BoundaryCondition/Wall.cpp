#include "Wall.hpp"

Compressible Wall::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = w;

    Vars<3> normalVector = vector3toVars(f.normalVector);
    Vars<3> ghostVelocity = w.velocity() - 2*w.normalVelocity(normalVector)*normalVector;
    double density = w.density();

    out[Compressible::RHO_U] = density*ghostVelocity[0];
    out[Compressible::RHO_V] = density*ghostVelocity[1];
    out[Compressible::RHO_W] = density*ghostVelocity[2];

    return out;
}

void Wall::correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (auto & faceIndex : boundary.facesIndex)
    {
        Vars<5> wlDiff = dot(grad[ownerIndexList[faceIndex]], vector3toVars(faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center));

        wl[faceIndex] = w[ownerIndexList[faceIndex]] + phi[ownerIndexList[faceIndex]]*wlDiff;
        wl[faceIndex].setThermoVar(thermoModel->updateThermo(wl[faceIndex]));

        wr[faceIndex] = calculateState(wl[faceIndex], faces[faceIndex], thermoModel);
    }

    /*for (auto & faceIndex : boundary.facesIndex)
    {
        Vars<5> wlDiff = dot(grad[ownerIndexList[faceIndex]], vector3toVars(faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center));

        wl[faceIndex] = w[ownerIndexList[faceIndex]] + phi[ownerIndexList[faceIndex]]*wlDiff;
        wl[faceIndex].setThermoVar(thermoModel->updateThermo(wl[faceIndex]));

        Vars<5> wrDiff = dot(grad[ownerIndexList[faceIndex]], -vector3toVars(faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center));
        wr[faceIndex] = wr[faceIndex] + phi[ownerIndexList[faceIndex]]*wrDiff;
    }*/
}