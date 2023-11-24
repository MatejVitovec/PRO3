#include <set>
#include <iostream>
#include "FVMScheme.hpp"

#include "BoundaryCondition/Periodicity.hpp"

#include "outputCFD.hpp"

void FVMScheme::setCfl(double cfl_)
{
    cfl = cfl_;
}

void FVMScheme::setMaxIter(int maxIter_)
{
    maxIter = maxIter_;
}

void FVMScheme::setTargetError(double targetError_)
{
    targetError = targetError_;
}

void FVMScheme::setLocalTimeStep(bool localTimeStep_)
{
    localTimeStep = localTimeStep_;
}

double FVMScheme::getCfl() const
{
    return cfl;
}

int FVMScheme::getMaxIter() const
{
    return maxIter;
}

double FVMScheme::getTargetError() const
{
    return targetError;
}

bool FVMScheme::getTimeStepsettings() const
{
    return localTimeStep;
}

const Mesh& FVMScheme::getMesh() const
{
    return mesh;
}



void FVMScheme::setInitialConditions(Compressible initialCondition)
{
    w = Field<Compressible>(mesh.getCellsSize());
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = initialCondition;
    }    
}

void FVMScheme::setInitialConditionsPrimitive(Vars<5> initialCondition)
{
    Compressible CompressibleIC = thermo->primitiveToConservative(initialCondition);

    w = Field<Compressible>(mesh.getCellsSize());
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = CompressibleIC;
    }  
}

void FVMScheme::init()
{
    gradientScheme->init(mesh, boundaryConditionList);

    wl = Field<Compressible>(mesh.getFacesSize());
    wr = Field<Compressible>(mesh.getFacesSize());
    wlOld = Field<Compressible>(mesh.getFacesSize());
    wrOld = Field<Compressible>(mesh.getFacesSize());

    timeSteps = Field<double>(mesh.getCellsSize());
}

void FVMScheme::applyBoundaryConditions()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<Face>& faceList = mesh.getFaceList();

    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->apply(ownerIndexList, faceList, w, wr, wrOld, thermo.get());
    }
}

void FVMScheme::setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions)
{
    boundaryConditionList = std::move(boundaryConditions);
}

void FVMScheme::calculateWlWr()
{
    //Without reconstruction

    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < mesh.getFacesSize(); i++)
    {        
        wl[i] = w[ownerIndexList[i]];

        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            wr[i] = w[neighbour];
        }
    }    
}

void FVMScheme::initWrWlOld()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < mesh.getFacesSize(); i++)
    {        
        wlOld[i] = w[ownerIndexList[i]];
        wrOld[i] = w[ownerIndexList[i]];
    }    

    wlOld = thermo->updateField(wlOld, wlOld);
    wrOld = thermo->updateField(wrOld, wrOld);
}

void FVMScheme::reconstruct()
{
    Field<Mat<5,3>> grad = gradientScheme->calculateGradient(wl, wr, mesh);

    Field<Vars<5>> phi = limiter->calculateLimiter(wl, wr, grad, mesh);

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            Vars<5> wlDiff = dot(grad[ownerIndexList[i]], vector3toVars(faces[i].midpoint - cells[ownerIndexList[i]].center));
            Vars<5> wrDiff = dot(grad[neighborIndexList[i]], vector3toVars(faces[i].midpoint - cells[neighborIndexList[i]].center));

            wl[i] = w[ownerIndexList[i]] + phi[ownerIndexList[i]]*wlDiff;
            wr[i] = w[neighborIndexList[i]] + phi[neighborIndexList[i]]*wrDiff;
        }
    }

    //limit na minimalni hodnoty u steny
    /*for (auto & boundaryCondition : boundaryConditionList)
    {
        if (boundaryCondition->getType() == BoundaryCondition::WALL)
        {
            std::vector<int> facesIndex = boundaryCondition->getBoundary().facesIndex;

            for (int i = 0; i < facesIndex.size(); i++)
            {
                if (wl[facesIndex[i]][Compressible::RHO] < 0.1)
                {
                    
                }
            }
        }
    }*/

    //korekce 2radu u steny
    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->correct(w, wl, wr, wrOld, grad, phi, mesh, thermo.get());
    }

    wr = thermo->updateInetrnalFieldFaces(wr, wrOld, mesh);

    for (auto & boundaryCondition : boundaryConditionList)
    {
        if (boundaryCondition->getType() == BoundaryCondition::PERIODICITY)
        {
            std::vector<int> boundaryFaces = boundaryCondition->getBoundary().facesIndex;
            std::vector<int> associatedBoundaryFaces = static_cast<Periodicity*>(boundaryCondition.get())->getPeriodicityFacesIndex();

            Vector3 shift = static_cast<Periodicity*>(boundaryCondition.get())->getFaceShift();

            for (int i = 0; i < boundaryFaces.size(); i++)
            {
                Vars<5> wlDiff = dot(grad[ownerIndexList[boundaryFaces[i]]], vector3toVars(faces[boundaryFaces[i]].midpoint - cells[ownerIndexList[boundaryFaces[i]]].center));
                
                Vars<5> wrDiff = dot(grad[ownerIndexList[associatedBoundaryFaces[i]]], vector3toVars(faces[boundaryFaces[i]].midpoint - cells[ownerIndexList[associatedBoundaryFaces[i]]].center + shift));

                wl[boundaryFaces[i]] = w[ownerIndexList[boundaryFaces[i]]] + phi[ownerIndexList[boundaryFaces[i]]]*wlDiff;                
                wr[boundaryFaces[i]] = w[ownerIndexList[associatedBoundaryFaces[i]]] + phi[ownerIndexList[associatedBoundaryFaces[i]]]*wrDiff;

                wr[boundaryFaces[i]].setThermoVar(thermo->updateThermo(wr[boundaryFaces[i]], wrOld[boundaryFaces[i]]));
            }
        }
    }
    
    wl = thermo->updateField(wl, wlOld);

    wlOld = wl;
    wrOld = wr;
}

void FVMScheme::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    for (int i = 0; i < w.size(); i++)
    {
        Vars<3> projectedArea = vector3toVars(cells[i].projectedArea);
        timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(w[i].velocity() + Vars<3>(w[i].soundSpeed()))));
    }

    if (!localTimeStep)
    {
        double timeStep = min(timeSteps);
        timeSteps = Field<double>(timeSteps.size(), timeStep);

        time += timeStep;
    }
}

void FVMScheme::calculateFluxes()
{
    fluxes = fluxSolver->claculateFluxes(wl, wr, mesh.getFaceList());
}

Field<Compressible> FVMScheme::getResults() const
{
    return w;
}