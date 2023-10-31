#include "BarthJespersen.hpp"


Field<Vars<5>> BarthJespersen::calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<std::array<Vars<3>, 5>>& grad, const Mesh& mesh) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Vars<5>> out(cells.size());

    Field<Vars<5>> wCmax(cells.size());
    Field<Vars<5>> wCmin(cells.size());
    //Field<Vars<5>> wCmax(cells.size(), std::array<double, 5>({0.0, 0.0, 0.0, 0.0, 0.0}));
    //Field<Vars<5>> wCmin(cells.size(), std::array<double, 5>({0.0, 0.0, 0.0, 0.0, 0.0}));

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

        for (int j = 0; j < cells[i].ownFaceIndex.size(); j++)
        {
            Vars<5> phiCn;
            Vars<5> denominator;            

            Vars<3> cellToFaceDist = vector3toVars(faces[cells[i].ownFaceIndex[j]].midpoint - cells[i].center);

            for (int k = 0; k < 5; k++)
            {
                denominator[k] = dot(grad[i][k], cellToFaceDist);
            }

            for (int k = 0; k < 5; k++)
            {
                if (denominator[k] > 0.0)
                {
                    phiCn[k] = std::min(1.0, std::max(0.0, (wCmax[i][k] - wC[k])/denominator[k]));
                }
                else if (denominator[k] < 0.0)
                {
                    phiCn[k] = std::min(1.0, std::max(0.0, (wCmin[i][k] - wC[k])/denominator[k]));
                }
                else
                {
                    phiCn[k] = 1.0;
                }                
            }

            phiC = min(phiC, phiCn);
        }

        for (int j = 0; j < cells[i].neighborFaceIndex.size(); j++)
        {
            Vars<5> phiCn;
            Vars<5> denominator;

            Vars<3> cellToFaceDist = vector3toVars(faces[cells[i].neighborFaceIndex[j]].midpoint - cells[i].center);

            for (int k = 0; k < 5; k++)
            {
                denominator[k] = dot(grad[i][k], cellToFaceDist);
            }

            for (int k = 0; k < 5; k++)
            {
                if (denominator[k] > 0.0)
                {
                    phiCn[k] = std::min(1.0, std::max(0.0, (wCmax[i][k] - wC[k])/denominator[k]));
                }
                else if (denominator[k] < 0.0)
                {
                    phiCn[k] = std::min(1.0, std::max(0.0, (wCmin[i][k] - wC[k])/denominator[k]));
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
}