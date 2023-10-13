#include <iostream>
#include "Periodicity.hpp"


void Periodicity::init(const Mesh& mesh)
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    
    periodicityFacesIndex.clear();
    periodicityFacesOwnersIndexes.clear();

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        Vector3 associatedFaceMidpoint = faceList[boundary.facesIndex[i]].midpoint;

        int j;
        bool isFound = false;
        for (j = 0; j < boundary.facesIndex.size(); j++)
        {
            if(associatedFaceMidpoint == (faceList[boundary.facesIndex[j]].midpoint + faceMidpointShift))
            {
                isFound = true;
                break;
            }
        }
        
        if(isFound)
        {
            periodicityFacesIndex.push_back(boundary.facesIndex[j]);
            periodicityFacesOwnersIndexes.push_back(ownerIndexList[boundary.facesIndex[j]]);
        }
        else
        {
            std::cout << "Periodicity face not found" << std::endl;
        }        
    }    
}


Compressible Periodicity::calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const
{
    //je nutne defonovat z duvodu abstraktni virtualni funkce - redefinuji velou funkci apply
    std::cout << "ERROR" <<std::endl; 
    return Compressible();
}

void Periodicity::apply(const std::vector<int>& ownerIndexList,const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const
{
    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        wr[boundary.facesIndex[i]] = w[ownerIndexList[i]];
    }   
}