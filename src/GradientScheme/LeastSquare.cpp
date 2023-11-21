#include "LeastSquare.hpp"


void LeastSquare::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<Mat<3,3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            for (int k = 0; k < 3; k++)
            {
                for (int kk = 0; kk < 3; kk++)
                {
                    M[i][k][kk] += cellToCellDelta[idx][kk]*cellToCellDelta[idx][k];
                }                
            }
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            for (int k = 0; k < 3; k++)
            {
                for (int kk = 0; kk < 3; kk++)
                {
                    M[i][k][kk] += (-cellToCellDelta[idx][kk])*(-cellToCellDelta[idx][k]);
                }                
            }
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

    for (int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            b[i] += outerProd(cellToCellDelta[idx], wr[idx] - wl[idx]);
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            b[i] += outerProd(-cellToCellDelta[idx], wl[idx] - wr[idx]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}