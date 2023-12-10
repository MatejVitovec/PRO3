#include "LeastSquare.hpp"


void LeastSquare::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<Mat<3,3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    /*const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();*/

    for (int i = 0; i < cells.size(); i++)
    {
        std::vector<int> cellFacesIndexes = cells[i].ownFaceIndex;
        cellFacesIndexes.insert(cellFacesIndexes.end(), cells[i].neighborFaceIndex.begin(), cells[i].neighborFaceIndex.end());

        for (int j = 0; j < cellFacesIndexes.size(); j++)
        {
            M[i] += outerProd(cellToCellDelta[cellFacesIndexes[j]], cellToCellDelta[cellFacesIndexes[j]]);
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
        std::vector<int> cellFacesIndexes = cells[i].ownFaceIndex;
        cellFacesIndexes.insert(cellFacesIndexes.end(), cells[i].neighborFaceIndex.begin(), cells[i].neighborFaceIndex.end());

        for (int j = 0; j < cellFacesIndexes.size(); j++)
        {
            b[i] += outerProd(cellToCellDelta[cellFacesIndexes[j]], wr[cellFacesIndexes[j]] - wl[cellFacesIndexes[j]]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}

Field<Mat<5,3>> LeastSquare::calculateGradient(const Field<Primitive>& ul, const Field<Primitive>& ur, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(ul.size());
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Mat<3,5>> b = Field<Mat<3,5>>(MInv.size());

    for (int i = 0; i < cells.size(); i++)
    {
        std::vector<int> cellFacesIndexes = cells[i].ownFaceIndex;
        cellFacesIndexes.insert(cellFacesIndexes.end(), cells[i].neighborFaceIndex.begin(), cells[i].neighborFaceIndex.end());

        for (int j = 0; j < cellFacesIndexes.size(); j++)
        {
            b[i] += outerProd(cellToCellDelta[cellFacesIndexes[j]], ur[cellFacesIndexes[j]] - ul[cellFacesIndexes[j]]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}