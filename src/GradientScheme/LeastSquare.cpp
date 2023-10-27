#include "LeastSquare.hpp"


void LeastSquare::init(const Mesh& mesh)
{
    Field<std::array<Vars<3>, 3>> M(mesh.getCellsSize());

    Field<Vars<3>> delta = Field<Vars<3>>(mesh.getFacesSize());

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
                    M[i][k][kk] += delta[idx][k]*delta[idx][kk];
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
                    M[i][k][kk] -= delta[idx][k]*delta[idx][kk];
                }                
            }
        }
    }

    calculateInverseM(M);
}


void LeastSquare::calculatesDeltas(const Mesh& mesh)
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < faces.size(); i++)
    {
        if (neighbours[i] == -1)
        {
            delta[i][0] = 2*(faces[i].midpoint.x - cells[owners[i]].center.x);
            delta[i][1] = 2*(faces[i].midpoint.y - cells[owners[i]].center.y);
            delta[i][2] = 2*(faces[i].midpoint.z - cells[owners[i]].center.z);
        }
        else
        {
            delta[i][0] = cells[neighbours[i]].center.x - cells[owners[i]].center.x;
            delta[i][1] = cells[neighbours[i]].center.y - cells[owners[i]].center.y;
            delta[i][2] = cells[neighbours[i]].center.z - cells[owners[i]].center.z;
        }
    }
}


void LeastSquare::calculateInverseM(Field<std::array<Vars<3>, 3>> M)
{
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


Field<std::array<Vars<5>, 3>> LeastSquare::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<std::array<Vars<5>, 3>> grad(wl.size());
    

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<std::array<Vars<3>, 5>> b = Field<std::array<Vars<3>, 5>>();


    for (int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            for (int k = 0; k < 5; k++)
            {
                b[i][k] += (wr[idx][k] - wl[idx][k])*delta[idx];
            }
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            for (int k = 0; k < 5; k++)
            {
                b[i][k] += (wl[idx][k] - wr[idx][k])*delta[idx];
            }
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            Vars<5> aux;
            for (int j = 0; j < 5; j++)
            {            
                aux[j] = dot(MInv[i][k], b[i][j]);
            }
            grad[i][k] = aux;
        }     
    }

    return grad;
}