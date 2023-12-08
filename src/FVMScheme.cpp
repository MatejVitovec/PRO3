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

    timeSteps = Field<double>(mesh.getCellsSize());
}

void FVMScheme::applyBoundaryConditions()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<Face>& faceList = mesh.getFaceList();

    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->apply(ownerIndexList, faceList, w, wr, thermo.get());
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

void FVMScheme::reconstruct()
{
    Field<Mat<5,3>> grad = gradientScheme->calculateGradient(wl, wr, mesh);

    Field<Vars<5>> phi = limiter->calculateLimiter(wl, wr, grad, mesh);

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    //EOS INTERVAL PRESERVING
    /*for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < 0.11 || w[i].internalEnergy() < 2010000.0)
        {
            phi[i] = Vars<5>(0.0);
            std::cout << "Err, EOS range; rho: " << w[i].density() << " e: " << w[i].internalEnergy() << std::endl;
        }
    }*/
    

    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            Vars<5> wlDiff = dot(grad[ownerIndexList[i]], vector3toVars(faces[i].midpoint - cells[ownerIndexList[i]].center));
            Vars<5> wrDiff = dot(grad[neighborIndexList[i]], vector3toVars(faces[i].midpoint - cells[neighborIndexList[i]].center));

            wl[i] = w[ownerIndexList[i]] + phi[ownerIndexList[i]]*wlDiff;
            wr[i] = w[neighborIndexList[i]] + phi[neighborIndexList[i]]*wrDiff;

            //EOS INTERVAL PRESERVING
            /*if (wl[i][Compressible::RHO] < 0.11 || wr[i][Compressible::RHO] < 0.11 || wl[i][Compressible::RHO_E] < 211000.0 || wr[i][Compressible::RHO_E] < 211000.0)
            {
                wl[i] = w[ownerIndexList[i]];
                wr[i] = w[neighborIndexList[i]];
            }*/
        }
    }

    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->correct(w, wl, wr, grad, phi, mesh, thermo.get());
    }
    
    wl = thermo->updateField(wl);
    wr = thermo->updateInetrnalFieldFaces(wr, mesh);
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