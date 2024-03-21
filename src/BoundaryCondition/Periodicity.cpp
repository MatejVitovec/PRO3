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

    /*double minFaceSize = 100000.0;
    for (int i = 0; i < faceList.size(); i++)
    {
        if (minFaceSize > faceList[i].area)
        {
            minFaceSize = faceList[i].area;
        }        
    }

    double numTol = minFaceSize/2.0;*/
    double numTol = 0.000001;
    
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

Compressible Periodicity::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    //je nutne definovat z duvodu abstraktni virtualni funkce - redefinuji celou funkci apply
    std::cout << "ERROR" <<std::endl; 
    return Compressible();
}

void Periodicity::apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const
{
    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        //periodicity nema update thermo
        wr[boundary.facesIndex[i]] = w[periodicityFacesOwnersIndexes[i]];
    }   
}


void Periodicity::correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        Vars<5> wlDiff = dot(grad[ownerIndexList[boundary.facesIndex[i]]], vector3toVars(faces[boundary.facesIndex[i]].midpoint - cells[ownerIndexList[boundary.facesIndex[i]]].center));
                
        Vars<5> wrDiff = dot(grad[ownerIndexList[periodicityFacesIndex[i]]], vector3toVars(faces[boundary.facesIndex[i]].midpoint - cells[ownerIndexList[periodicityFacesIndex[i]]].center + faceMidpointShift));

        wl[boundary.facesIndex[i]] = w[ownerIndexList[boundary.facesIndex[i]]] + phi[ownerIndexList[boundary.facesIndex[i]]]*wlDiff;                
        wr[boundary.facesIndex[i]] = w[ownerIndexList[periodicityFacesIndex[i]]] + phi[ownerIndexList[periodicityFacesIndex[i]]]*wrDiff;

        wr[boundary.facesIndex[i]].setThermoVar(thermoModel->updateThermo(wr[boundary.facesIndex[i]]));
    }
}