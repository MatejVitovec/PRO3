#include <cmath>

#include "GradientScheme.hpp"

void GradientScheme::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{

}

Field<Mat<5,3>> GradientScheme::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(wl.size());

    return grad;
}

void GradientScheme::calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    cellToCellDelta = Field<Vars<3>>(faces.size());

    for (int i = 0; i < faces.size(); i++)
    {
        if (neighbours[i] != -1)
        {
            cellToCellDelta[i] = vector3toVars(cells[neighbours[i]].center - cells[owners[i]].center);
        }
    }

    for (auto & boundaryCondition : boundaryConditionList)
    {
        if (boundaryCondition->getType() == BoundaryCondition::PERIODICITY)
        {
            std::vector<int> boundaryFaces = boundaryCondition->getBoundary().facesIndex;
            std::vector<int> associatedBoundaryFaces = static_cast<Periodicity*>(boundaryCondition.get())->getPeriodicityFacesIndex();

            Vector3 faceMidpointShift = static_cast<Periodicity*>(boundaryCondition.get())->getFaceShift();

            for (int i = 0; i < boundaryFaces.size(); i++)
            {
                cellToCellDelta[boundaryFaces[i]] = vector3toVars(cells[owners[associatedBoundaryFaces[i]]].center - cells[owners[boundaryFaces[i]]].center - faceMidpointShift);
            }
        }
        else
        {
            std::vector<int> boundaryFaces = boundaryCondition->getBoundary().facesIndex;

            for (int i = 0; i < boundaryFaces.size(); i++)
            {
                //TODO - mozna to bude fungovat
                cellToCellDelta[boundaryFaces[i]] = vector3toVars(2*(faces[boundaryFaces[i]].midpoint - cells[owners[boundaryFaces[i]]].center));
            }
        }
    }
}
