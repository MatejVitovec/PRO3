#include "ExplicitEuler.hpp"
#include <iostream>
#include "outputCFD.hpp"


void ExplicitEuler::solve()
{
    init();

    w = thermo->updateField(w);

    int iter = 0;

    bool exitLoop = false;

    Field<Compressible> wn = Field<Compressible>(w.size());
    Vars<5> resNorm;

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        applyBoundaryConditions();
        //applyBoundaryConditionsPrimitive();

        calculateWlWr();
        //calculateUlUr();

        if(reconstruction)
        {
            reconstruct();
            //reconstructPrimitive();
        }

        calculateFluxes();

        Field<Vars<5>> res = calculateResidual();

        wn = w + (res*timeSteps);

        wn = thermo->updateField(wn);

        resNorm = (wn - w).norm();
        outputCFD::saveResidual("../results/residuals.txt", resNorm);

        if(resNorm[0] < targetError) exitLoop = true;
        
        w = wn;

        if(iter % 200 == 0)
        {
            outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
        }
    }

    outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
    std::cout << "iter: " << iter << std::endl;

    std::cout << "time: " << time << std::endl;
}