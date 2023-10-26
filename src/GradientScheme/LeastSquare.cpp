#include "LeastSquare.hpp"


void LeastSquare::init(const Mesh& mesh)
{
    Field<Vars<3>> xRowOfM = Field<Vars<3>>(mesh.getCellsSize());
    Field<Vars<3>> yRowOfM = Field<Vars<3>>(mesh.getCellsSize());
    Field<Vars<3>> zRowOfM = Field<Vars<3>>(mesh.getCellsSize());

    Field<Vars<3>> delta = Field<Vars<3>>(mesh.getFacesSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < cells.size(); i++)
    {
        xRowOfM[i][0] = 0;
        xRowOfM[i][1] = 0;
        xRowOfM[i][2] = 0;
        yRowOfM[i][0] = 0;
        yRowOfM[i][1] = 0;
        yRowOfM[i][2] = 0;
        zRowOfM[i][0] = 0;
        zRowOfM[i][1] = 0;
        zRowOfM[i][2] = 0;

        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            int idx = cells[i].ownFaceIndex[j];

            xRowOfM[i][0] += delta[idx][0]*delta[idx][0];
            xRowOfM[i][1] += delta[idx][0]*delta[idx][1];
            xRowOfM[i][2] += delta[idx][0]*delta[idx][2];

            yRowOfM[i][0] += delta[idx][0]*delta[idx][1];
            yRowOfM[i][1] += delta[idx][1]*delta[idx][1];
            yRowOfM[i][2] += delta[idx][1]*delta[idx][2];

            zRowOfM[i][0] += delta[idx][0]*delta[idx][2];
            zRowOfM[i][1] += delta[idx][1]*delta[idx][2];
            zRowOfM[i][2] += delta[idx][2]*delta[idx][2];
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            xRowOfM[i][0] += (-delta[idx][0])*(-delta[idx][0]);
            xRowOfM[i][1] += (-delta[idx][0])*(-delta[idx][1]);
            xRowOfM[i][2] += (-delta[idx][0])*(-delta[idx][2]);

            yRowOfM[i][0] += (-delta[idx][0])*(-delta[idx][1]);
            yRowOfM[i][1] += (-delta[idx][1])*(-delta[idx][1]);
            yRowOfM[i][2] += (-delta[idx][1])*(-delta[idx][2]);

            zRowOfM[i][0] += (-delta[idx][0])*(-delta[idx][2]);
            zRowOfM[i][1] += (-delta[idx][1])*(-delta[idx][2]);
            zRowOfM[i][2] += (-delta[idx][2])*(-delta[idx][2]);
        }
    }

    calculateInverseM(xRowOfM, yRowOfM, zRowOfM);

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


void LeastSquare::calculateInverseM(Field<Vars<3>> xRowOfM, Field<Vars<3>> yRowOfM, Field<Vars<3>> zRowOfM)
{
    for (int i = 0; i < xRowOfMInv.size(); i++)
    {
        double det = det3by3(xRowOfM[i], yRowOfM[i], zRowOfM[i]);
        std::array<double, 9> adj = adjoint3by3(xRowOfM[i], yRowOfM[i], zRowOfM[i]);

        xRowOfMInv[i][0] = adj[0]/det;
        xRowOfMInv[i][1] = -adj[1]/det;
        xRowOfMInv[i][2] = adj[2]/det;
        yRowOfMInv[i][0] = -adj[3]/det;
        yRowOfMInv[i][1] = adj[4]/det;
        yRowOfMInv[i][2] = -adj[5]/det;
        zRowOfMInv[i][0] = adj[6]/det;
        zRowOfMInv[i][1] = -adj[7]/det;
        zRowOfMInv[i][2] = adj[8]/det;
    }
}


double LeastSquare::det3by3(Vars<3> xRow, Vars<3> yRow, Vars<3> zRow) const
{
    return xRow[0]*yRow[1]*zRow[2] + xRow[1]*yRow[2]*zRow[0] + xRow[2]*yRow[0]*zRow[1] - xRow[0]*yRow[2]*zRow[1] - xRow[1]*yRow[0]*zRow[2] - xRow[2]*yRow[1]*zRow[0];
}


std::array<double, 9> LeastSquare::adjoint3by3(Vars<3> xRow, Vars<3> yRow, Vars<3> zRow) const
{
    return {yRow[1]*zRow[2] - yRow[2]*zRow[1],
            yRow[0]*zRow[2] - yRow[2]*zRow[0],
            yRow[0]*zRow[1] - yRow[1]*zRow[0],
            xRow[1]*zRow[2] - xRow[2]*zRow[1],
            xRow[0]*zRow[2] - xRow[2]*zRow[0],
            xRow[0]*zRow[1] - xRow[1]*zRow[0],
            xRow[1]*yRow[2] - xRow[2]*yRow[1],
            xRow[0]*yRow[2] - xRow[2]*yRow[0],
            xRow[0]*yRow[1] - xRow[1]*yRow[0]};
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

            b[i][0] += (wr[idx][0] - wl[idx][0])*delta[idx];
            b[i][1] += (wr[idx][1] - wl[idx][1])*delta[idx];
            b[i][2] += (wr[idx][2] - wl[idx][2])*delta[idx];
            b[i][3] += (wr[idx][3] - wl[idx][3])*delta[idx];
            b[i][4] += (wr[idx][4] - wl[idx][4])*delta[idx];
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            int idx = cells[i].neighborFaceIndex[j];

            b[i][0] += (wl[idx][0] - wr[idx][0])*delta[idx];
            b[i][1] += (wl[idx][1] - wr[idx][1])*delta[idx];
            b[i][2] += (wl[idx][2] - wr[idx][2])*delta[idx];
            b[i][3] += (wl[idx][3] - wr[idx][3])*delta[idx];
            b[i][4] += (wl[idx][4] - wr[idx][4])*delta[idx];
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        Vars<5> aux;
        for (int j = 0; j < 5; j++)
        {            
            aux[j] = dot(xRowOfMInv[i], b[i][j]);
        }
        grad[i][0] = aux;

        for (int j = 0; j < 5; j++)
        {            
            aux[j] = dot(yRowOfMInv[i], b[i][j]);
        }
        grad[i][1] = aux;

        for (int j = 0; j < 5; j++)
        {            
            aux[j] = dot(yRowOfMInv[i], b[i][j]);
        }
        grad[i][2] = aux;        
    }

    return grad;
}