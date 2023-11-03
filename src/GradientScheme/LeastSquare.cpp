#include "LeastSquare.hpp"


void LeastSquare::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<std::array<Vars<3>, 3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < cells.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            for (int kk = 0; kk < 3; kk++)
            {
                M[i][k][kk] = 0;
            }            
        }

        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            for (int k = 0; k < 3; k++)
            {
                for (int kk = 0; kk < 3; kk++)
                {
                    M[i][k][kk] += cellToCellDelta[idx][k]*cellToCellDelta[idx][kk];
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
                    M[i][k][kk] -= cellToCellDelta[idx][k]*cellToCellDelta[idx][kk];
                }                
            }
        }
    }

    calculateInverseM(M);
}



void LeastSquare::calculateInverseM(Field<std::array<Vars<3>, 3>> M)
{
    MInv = Field<std::array<Vars<3>, 3>>(M.size());

    for (int i = 0; i < MInv.size(); i++)
    {
        double det = det3by3(M[i]);
        std::array<Vars<3>, 3> adj = adjoint3by3(M[i]);

        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if(j*k % 2 == 0)
                    MInv[i][j][k] = adj[j][k]/det;
                else
                    MInv[i][j][k] = -adj[j][k]/det;                
            }            
        }
    }
}


double LeastSquare::det3by3(std::array<Vars<3>, 3> M) const
{
    return M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[0][0]*M[1][2]*M[2][1] - M[0][1]*M[1][0]*M[2][2] - M[0][2]*M[1][1]*M[2][0];
}


std::array<Vars<3>, 3> LeastSquare::adjoint3by3(std::array<Vars<3>, 3> M) const
{
    return {Vars<3>({M[1][1]*M[2][2] - M[1][2]*M[2][1], M[1][0]*M[2][2] - M[1][2]*M[2][0], M[1][0]*M[2][1] - M[1][1]*M[2][0]}),
            Vars<3>({M[0][1]*M[2][2] - M[0][2]*M[2][1], M[0][0]*M[2][2] - M[0][2]*M[2][0], M[0][0]*M[2][1] - M[0][1]*M[2][0]}),
            Vars<3>({M[0][1]*M[1][2] - M[0][2]*M[1][1], M[0][0]*M[1][2] - M[0][2]*M[1][0], M[0][0]*M[1][1] - M[0][1]*M[1][0]})};
}


Field<std::array<Vars<3>, 5>> LeastSquare::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<std::array<Vars<3>, 5>> grad(wl.size());
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<std::array<Vars<3>, 5>> b = Field<std::array<Vars<3>, 5>>(MInv.size());

    for (int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            Vars<5> aux = wr[idx] - wl[idx];

            for (int k = 0; k < 5; k++)
            {
                b[i][k] += aux[k]*cellToCellDelta[idx];
            }
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            Vars<5> aux = wl[idx] - wr[idx];

            for (int k = 0; k < 5; k++)
            {
                b[i][k] -= aux[k]*cellToCellDelta[idx];
            }
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {

        for (int j = 0; j < 5; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                 grad[i][j][k] = dot(MInv[i][k], b[i][j]);
            }
        }
        
    }

    return grad;
}