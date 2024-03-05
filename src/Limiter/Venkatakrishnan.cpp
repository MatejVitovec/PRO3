#include "Venkatakrishnan.hpp"
#include <iostream>

/*Field<Vars<5>> Venkatakrishnan::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Vars<5>> out(cells.size());

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());

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

    for (int i = 0; i < cells.size(); i++)
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

            Vars<5>denominator = dot(grad[i], cellToFaceDist);

            for (int k = 0; k < 5; k++)
            {
                double y;
                if (denominator[k] > 0.0)
                {
                    y = std::max(0.0, (wCmax[i][k] - wC[k])/denominator[k]);
                    phiCn[k] = (y*y + 2*y)/(y*y + y + 2);
                }
                else if (denominator[k] < 0.0)
                {
                    y = std::max(0.0, (wCmin[i][k] - wC[k])/denominator[k]);
                    phiCn[k] = (y*y + 2*y)/(y*y + y + 2);
                }
                else
                {
                    phiCn[k] = 1.0;
                }                
            }

            phiC = min(phiC, phiCn);
        }

        out[i] = phiC;
    }

    return out;
}*/

/*Field<Vars<5>> Venkatakrishnan::calculateLimiter(const Field<Primitive>& ul, const Field<Primitive>& ur, const Field<Mat<5,3>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Vars<5>> out(cells.size());

    Field<Vars<5>> uCmax(cells.size());
    Field<Vars<5>> uCmin(cells.size());

    for (int i = 0; i < faces.size(); i++)
    {
        int ownerIndex = owners[i];  

        uCmax[ownerIndex] = max(uCmax[ownerIndex], ur[i]);
        uCmin[ownerIndex] = min(uCmin[ownerIndex], ur[i]);

        int neighbourIndex = neighbours[i];
        if (neighbourIndex >= 0)
        {
            uCmax[neighbourIndex] = max(uCmax[neighbourIndex], ul[i]);
            uCmin[neighbourIndex] = min(uCmin[neighbourIndex], ul[i]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {        
        Vars<5> phiC({10.0, 10.0, 10.0, 10.0, 10.0});

        Primitive uC;
        if(cells[i].ownFaceIndex.size() > 0)
            uC = ul[cells[i].ownFaceIndex[0]];
        else
            uC = ur[cells[i].neighborFaceIndex[0]];

        std::vector<int> cellFacesIndexes = cells[i].ownFaceIndex;
        cellFacesIndexes.insert(cellFacesIndexes.end(), cells[i].neighborFaceIndex.begin(), cells[i].neighborFaceIndex.end());

        for (int j = 0; j < cellFacesIndexes.size(); j++)
        {
            Vars<5> phiCn;

            Vars<3> cellToFaceDist = vector3toVars(faces[cellFacesIndexes[j]].midpoint - cells[i].center);

            Vars<5>denominator = dot(grad[i], cellToFaceDist);

            for (int k = 0; k < 5; k++)
            {
                double y;
                if (denominator[k] > 0.0)
                {
                    y = std::max(0.0, (uCmax[i][k] - uC[k])/denominator[k]);
                    phiCn[k] = (y*y + 2*y)/(y*y + y + 2);
                }
                else if (denominator[k] < 0.0)
                {
                    y = std::max(0.0, (uCmin[i][k] - uC[k])/denominator[k]);
                    phiCn[k] = (y*y + 2*y)/(y*y + y + 2);
                }
                else
                {
                    phiCn[k] = 1.0;
                }                
            }

            phiC = min(phiC, phiCn);
        }

        out[i] = phiC;
    }

    return out;
}*/

double Venkatakrishnan::limiterFunction(double y) const
{
    return (y*y + 2*y)/(y*y + y + 2);
}