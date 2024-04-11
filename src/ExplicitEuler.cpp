#include "ExplicitEuler.hpp"
#include <iostream>
#include "outputCFD.hpp"


void ExplicitEuler::solve()
{
    init();

    double lastTmdUpdateTime = 0.0;
    long thermoUpdateFieldTime = 0.0;

    w = thermo->updateField(w);

    int iter = 0;

    bool exitLoop = false;

    Field<Compressible> wn = Field<Compressible>(w.size());
    Vars<5> resNorm;

    auto startAll = std::chrono::high_resolution_clock::now();

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        applyBoundaryConditions();

        calculateWlWr();

        if(reconstruction)
        {
            reconstruct();
        }

        calculateFluxes();

        Field<Vars<5>> res = calculateResidual();

        wn = w + (res*timeSteps);

        auto start = std::chrono::high_resolution_clock::now();

        wn = thermo->updateField(wn);
        
        auto stop = std::chrono::high_resolution_clock::now();
        thermoUpdateFieldTime += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

        if(iter % 100 == 0)
        {
            outputCFD::saveValue("../results/time.txt", std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startAll).count());

            resNorm = res.norm();
            outputCFD::saveResidual("../results/residuals.txt", resNorm);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;

            outputCFD::saveValue("../results/tmdUpdateTime.txt", thermoUpdateFieldTime/1000.0 - lastTmdUpdateTime);
            lastTmdUpdateTime = thermoUpdateFieldTime/1000.0;
        }

        w = wn;

        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
        }
    }

    outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
    std::cout << "iter: " << iter << std::endl;

    std::cout << "time: " << time << std::endl;

    std::cout << "timeForCalculation: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startAll).count() << std::endl;
}