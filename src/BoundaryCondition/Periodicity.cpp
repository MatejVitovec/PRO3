#include <iostream>
#include "Periodicity.hpp"


Periodicity::Periodicity(Boundary meshBoundary, Vector3 faceMidpointShift_, std::string associatedBoundaryName_, const Mesh& mesh) : BoundaryCondition(meshBoundary, PERIODICITY), faceMidpointShift(faceMidpointShift_), associatedBoundaryName(associatedBoundaryName_)
{
    init(mesh);
}

void Periodicity::init(const Mesh& mesh)
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    double minFaceSize = 100000.0;
    for (int i = 0; i < faceList.size(); i++)
    {
        if (minFaceSize > faceList[i].area)
        {
            minFaceSize = faceList[i].area;
        }        
    }

    double numTol = minFaceSize/2.0;
    
    
    periodicityFacesIndex.clear();
    periodicityFacesOwnersIndexes.clear();

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        Vector3 associatedFaceMidpoint = faceList[boundary.facesIndex[i]].midpoint + faceMidpointShift;

        int j;
        bool isFound = false;
        for (j = 0; j < faceList.size(); j++)
        {
            if(norm2(associatedFaceMidpoint - faceList[j].midpoint) < numTol)
            {
                isFound = true;
                break;
            }
        }
        
        if(isFound)
        {
            periodicityFacesIndex.push_back(j);
            periodicityFacesOwnersIndexes.push_back(ownerIndexList[j]);
        }
        else
        {
            std::cout << "Periodicity face not found" << std::endl;
        }        
    }    
}

std::vector<int> Periodicity::getPeriodicityFacesIndex() const
{
    return periodicityFacesIndex;
}

Vector3 Periodicity::getFaceShift() const
{
    return faceMidpointShift;
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
        //wr[boundary.facesIndex[i]] = w[ownerIndexList[i]];
        wr[boundary.facesIndex[i]] = w[periodicityFacesOwnersIndexes[i]];        
    }   
}