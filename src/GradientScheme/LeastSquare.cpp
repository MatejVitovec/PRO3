#include "LeastSquare.hpp"
#include <iostream>

void LeastSquare::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<Mat<3,3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < owners.size(); i++)
    {
        Mat<3,3> auxM = outerProd(cellToCellDelta[i], cellToCellDelta[i]);
        M[owners[i]] += auxM;
        if (neighbours[i] >= 0)
        {
            M[neighbours[i]] += auxM;
        }
    }

    calculateInverseM(M);
}


void LeastSquare::calculateInverseM(Field<Mat<3,3>> M)
{
    MInv = Field<Mat<3,3>>(M.size());

    for (int i = 0; i < MInv.size(); i++)
    {
        MInv[i] = inv(M[i]);
    }
}


Field<Mat<5,3>> LeastSquare::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(wl.size());
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Mat<3,5>> b = Field<Mat<3,5>>(MInv.size());

    for (int i = 0; i < owners.size(); i++)
    {
        int neighbour = neighbours[i];

        Mat<3,5> rhs = outerProd(cellToCellDelta[i], wr[i] - wl[i]);
        b[owners[i]] += rhs;
        if (neighbour >= 0)
        {
            b[neighbour] -= rhs;
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}