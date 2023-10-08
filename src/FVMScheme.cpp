#include "FVMScheme.hpp"


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

void FVMScheme::setInitialConditionsRiemann(Compressible initialConditionL, Compressible initialConditionR)
{
    w = Field<Compressible>(mesh.getCellsSize());
    int i = 0;
    for (i; i < 30; i++)
    {
        w[i] = initialConditionL;
    }
    for (i; i < mesh.getCellsSize(); i++)
    {
        w[i] = initialConditionR;
    }
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

void FVMScheme::setBoundaryConditions(std::vector<std::unique_ptr<BoundaryCondition>> boundaryConditions)
{
    boundaryConditionList = std::move(boundaryConditions);
}

void FVMScheme::calculateWlWr()
{
    //Without reconstruction

    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < faces.size(); i++)
    {        
        wl[i] = w[ownerIndexList[i]];

        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            wr[i] = w[neighbour];
        }
    }    
}

void FVMScheme::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    timeStep = 1000000;

    for (int i = 0; i < w.size(); i++)
    {
        double soundSpeed = w[i].soundSpeed();
        
        timeStep = std::min(cfl*((cells[i].volume)/(cells[i].projectedArea.x*(w[i].velocityU() + soundSpeed) + cells[i].projectedArea.y*(w[i].velocityV() + soundSpeed) + cells[i].projectedArea.z*(w[i].velocityW() + soundSpeed))), timeStep);
    }
    time += timeStep;
}

void FVMScheme::calculateFluxes()
{
    fluxes = fluxSolver->claculateFluxes(wl, wr, mesh.getFaceList());
}

Field<Compressible> FVMScheme::getResults() const
{
    return w;
}