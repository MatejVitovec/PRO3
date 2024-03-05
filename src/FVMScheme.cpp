#include <set>
#include <iostream>
#include "FVMScheme.hpp"

#include "BoundaryCondition/Periodicity.hpp"

#include "outputCFD.hpp"

void FVMScheme::setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_)
{
    gradientScheme = std::move(gradScheme_);
}
void FVMScheme::setReconstructionLimiter(std::unique_ptr<Limiter> limiter_)
{
    limiter = std::move(limiter_);
}

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

void FVMScheme::setReconstructionSettings(bool reconstruction_)
{
    reconstruction = reconstruction_;
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

bool FVMScheme::getReconstructionSettings() const
{
    return reconstruction;
}

const Mesh& FVMScheme::getMesh() const
{
    return mesh;
}

const Thermo* FVMScheme::getThermoRef()
{
    return thermo.get();
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
    ul = Field<Primitive>(mesh.getFacesSize());
    ur = Field<Primitive>(mesh.getFacesSize());

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

void FVMScheme::applyBoundaryConditionsPrimitive()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<Face>& faceList = mesh.getFaceList();

    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->apply(ownerIndexList, faceList, w, ur, thermo.get());
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

void FVMScheme::calculateUlUr()
{
    //Without reconstruction

    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < mesh.getFacesSize(); i++)
    {        
        ul[i] = Primitive(w[ownerIndexList[i]]);

        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            ur[i] = Primitive(w[neighbour]);
        }
    }    
}

void FVMScheme::boundField()
{
    constexpr double minDensity = 0.1;
    constexpr double maxMach = 2.2;
    constexpr double minInternalEnergy = 2100000.0;
    constexpr double maxInternalEnergy = 2800000.0;
    

    const std::vector<Cell>& cells = mesh.getCellList();
    for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < minDensity)
        {
            //std::cout << "rho: " << w[i].density() << " U: " << w[i].absVelocity() << " energy: " << w[i].internalEnergy() << std::endl;
            w[i][Compressible::RHO_U] = w[i].velocityU()*minDensity;
            w[i][Compressible::RHO_V] = w[i].velocityV()*minDensity;
            w[i][Compressible::RHO_W] = w[i].velocityW()*minDensity;
            w[i][Compressible::RHO_E] = w[i].totalEnergy()*minDensity;
            w[i][Compressible::RHO] = minDensity;
        }

        /*if (w[i].machNumber() > maxMach)
        {
            double velocityCoeff = (maxMach*w[i].soundSpeed())/(w[i].absVelocity());
            w[i][Compressible::RHO_U] = w[i][Compressible::RHO_U]*velocityCoeff;
            w[i][Compressible::RHO_V] = w[i][Compressible::RHO_V]*velocityCoeff;
            w[i][Compressible::RHO_W] = w[i][Compressible::RHO_W]*velocityCoeff;
            //std::cout << "M i: " << std::endl;
        }*/

        if (w[i].internalEnergy() < minInternalEnergy)
        {
            w[i][Compressible::RHO_E] = minInternalEnergy*w[i].density();
            //std::cout << "e min i: " << i << std::endl;
        }

        /*if (w[i].internalEnergy() > maxInternalEnergy)
        {
            w[i][Compressible::RHO_E] = maxInternalEnergy*w[i].density();
            //std::cout << "e min i: " << i << std::endl;
        }*/
    }
}

void FVMScheme::reconstruct()
{
    Field<Mat<5,3>> grad = gradientScheme->calculateGradient(wl, wr, mesh);

    Field<Vars<5>> phi = limiter->calculateLimiter(wl, wr, grad, mesh);

    //phi = phi*0.8;

    //phi[386] = Vars<5>(0.0);
    //phi[7385] = Vars<5>(0.0);

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    //EOS INTERVAL PRESERVING
    for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < 0.11 || w[i].internalEnergy() < 2200000.0)
        {
            phi[i] = Vars<5>(0.0);
            //std::cout << "Err, EOS range; rho: " << w[i].density() << " e: " << w[i].internalEnergy() << std::endl;
        }
    }
    

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

void FVMScheme::reconstructPrimitive()
{
    Field<Mat<5,3>> grad = gradientScheme->calculateGradient(ul, ur, mesh);

    Field<Vars<5>> phi = limiter->calculateLimiter(ul, ur, grad, mesh);

    //phi[386] = Vars<5>(0.0);
    //phi[7385] = Vars<5>(0.0);

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    /*for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < 0.1 || w[i].internalEnergy() < 2200000.0)
        {
            phi[i] = Vars<5>(0.0);
        }
    }*/

    for (int i = 0; i < faces.size(); i++)
    {
        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            Vars<5> ulDiff = dot(grad[ownerIndexList[i]], vector3toVars(faces[i].midpoint - cells[ownerIndexList[i]].center));
            Vars<5> urDiff = dot(grad[neighborIndexList[i]], vector3toVars(faces[i].midpoint - cells[neighborIndexList[i]].center));

            ul[i] = ul[i] + phi[ownerIndexList[i]]*ulDiff;
            ur[i] = ur[i] + phi[neighborIndexList[i]]*urDiff;
        }
    }

    for (auto & boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->correct(w, ul, ur, grad, phi, mesh, thermo.get());
    }
    
    ul = thermo->updateField(ul);
    ur = thermo->updateInetrnalFieldFaces(ur, mesh);
}

void FVMScheme::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    for (int i = 0; i < w.size(); i++)
    {
        Vars<3> projectedArea = vector3toVars(cells[i].projectedArea);
        timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(abs(w[i].velocity()) + Vars<3>(w[i].soundSpeed()))));
        if (timeSteps[i] < 0.0)
        {
            std::cout << "cell index: " << i << std::endl;
        }
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
    fluxes = fluxSolver->calculateFluxes(wl, wr, mesh.getFaceList());
    //fluxes = fluxSolver->calculateFluxes(ul, ur, mesh.getFaceList());
}

Field<Vars<5>> FVMScheme::calculateResidual()
{
    //TODO predelat na cyklus pres steny
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();

    Field<Vars<5>> res(w.size());

    for (int i = 0; i < cells.size(); i++)
    {
        Vars<5> aux= Vars<5>();

        for (auto & faceIndex : cells[i].ownFaceIndex)
        {
            aux -= fluxes[faceIndex];
        }

        for (auto & faceIndex : cells[i].neighborFaceIndex)
        {
            aux += fluxes[faceIndex];
        }

        res[i] = aux/(cells[i].volume);
    }    
    
    return res;
}

Field<Compressible> FVMScheme::getResults() const
{
    return w;
}