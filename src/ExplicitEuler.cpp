#include "ExplicitEuler.hpp"
#include <iostream>
#include "outputCFD.hpp"


void ExplicitEuler::solve()
{
    init();

    Field<Compressible> wOld = Field<Compressible>(w.size());

    w = thermo->updateField(w, w);
    initWrWlOld();

    int iter = 0;

    bool exitLoop = false;

    Field<Compressible> wn = Field<Compressible>(w.size());
    Vars<5> resNorm;

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        updateTimeStep();

        applyBoundaryConditions();

        calculateWlWr();

        reconstruct();

        calculateFluxes();

        Field<Vars<5>> res = calculateResidual();
        
        wn = w + (res*timeSteps);

        wn = thermo->updateField(wn, w);

        resNorm = (wn - w).norm();
        outputCFD::saveResidual("../results/residuals.txt", resNorm);

        if(resNorm[0] < targetError) exitLoop = true;
        
        w = wn;

        if(iter % 200 == 0 || iter > 8300)
        {
            outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
        }
    }

    outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
    std::cout << "iter: " << iter << std::endl;

    std::cout << "time: " << time << std::endl;
}

Field<Vars<5>> ExplicitEuler::calculateResidual()
{
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