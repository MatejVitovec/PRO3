#include <cmath>
#include <iostream>

#include "Limiter.hpp"

Field<Vars<5>> Limiter::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    //Field<Vars<5>> out(cells.size());
    Field<Vars<5>> out(cells.size(), Vars<5>({10.0, 10.0, 10.0, 10.0, 10.0}));

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());

    //over faces
    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];  

        wCmax[ownerIndex] = max(wCmax[ownerIndex], wr[i]);
        wCmin[ownerIndex] = min(wCmin[ownerIndex], wr[i]);

        int neighbourIndex = neighbours[i];
        if (neighbourIndex >= 0)
        {
            wCmax[neighbourIndex] = max(wCmax[neighbourIndex], wl[i]);
            wCmin[neighbourIndex] = min(wCmin[neighbourIndex], wl[i]);
        }
    }

    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbours[i];

        Compressible wCOwner = wl[i];
        Vars<5> denominatorOwner = dot(grad[owner], vector3toVars(faces[i].midpoint - cells[owner].center));
        Vars<5> phiCnOwner;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorOwner[k] > 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmax[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else if (denominatorOwner[k] < 0.0)
            {
                phiCnOwner[k] = limiterFunction(std::max(0.0, (wCmin[owner][k] - wCOwner[k])/denominatorOwner[k]));
            }
            else
            {
                phiCnOwner[k] = 1.0;
            }                
        }

        out[owner] = min(out[owner], phiCnOwner);

        if (neighbor < 0)
        {
            continue;
        }
        
        Compressible wCNeighbor = wr[i];
        Vars<5> denominatorNeighbor = dot(grad[neighbor], vector3toVars(faces[i].midpoint - cells[neighbor].center));
        Vars<5> phiCnNeighbor;

        for (int k = 0; k < 5; k++)
        {
            if (denominatorNeighbor[k] > 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmax[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else if (denominatorNeighbor[k] < 0.0)
            {
                phiCnNeighbor[k] = limiterFunction(std::max(0.0, (wCmin[neighbor][k] - wCNeighbor[k])/denominatorNeighbor[k]));
            }
            else
            {
                phiCnNeighbor[k] = 1.0;
            }                
        }

        out[neighbor] = min(out[neighbor], phiCnNeighbor);
    }

    //over cells
    /*for (int i = 0; i < cells.size(); i++)
    {        
        Vars<5> phiC({10.0, 10.0, 10.0, 10.0, 10.0});

        Compressible wC;
        if(cells[i].ownFaceIndex.size() > 0)
            wC = wl[cells[i].ownFaceIndex[0]];
        else
            wC = wr[cells[i].neighborFaceIndex[0]];


        std::vector<int> cellFacesIndexes = cells[i].ownFaceIndex;
        cellFacesIndexes.insert(cellFacesIndexes.end(), cells[i].neighborFaceIndex.begin(), cells[i].neighborFaceIndex.end());

        for (int j = 0; j < cellFacesIndexes.size(); j++)
        {
            Vars<5> phiCn;

            Vars<3> cellToFaceDist = vector3toVars(faces[cellFacesIndexes[j]].midpoint - cells[i].center);

            Vars<5> denominator = dot(grad[i], cellToFaceDist);

            for (int k = 0; k < 5; k++)
            {
                if (denominator[k] > 0.0)
                {
                    phiCn[k] = limiterFunction(std::max(0.0, (wCmax[i][k] - wC[k])/denominator[k]));
                }
                else if (denominator[k] < 0.0)
                {
                    phiCn[k] = limiterFunction(std::max(0.0, (wCmin[i][k] - wC[k])/denominator[k]));
                }
                else
                {
                    phiCn[k] = 1.0;
                }                
            }

            phiC = min(phiC, phiCn);
        }

        out[i] = phiC;
    }*/

    return out;
}

double Limiter::limiterFunction(double y) const
{
    return 0.0;
}