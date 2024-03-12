#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>

#include "outputCFD.hpp"
#include "CaseSetter.hpp"


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    CaseSetter setter = CaseSetter();
    setter.loadSettingFile("../case/setup.txt");

    std::unique_ptr<FVMScheme> mySolver = setter.createAndSetSolver();

    outputCFD::outputVTK("../results/results.0.vtk", mySolver->getMesh(), mySolver->getResults());

    auto stop1 = std::chrono::high_resolution_clock::now();

    mySolver->solve();
    
    auto stop2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1).count() << " ms\n";

    outputCFD::saveFieldOnBoundary("../results/pressure.txt", "wall", mySolver->getMesh(), mySolver->getResults());

    outputCFD::outputVTKPeriodicBoundary("../results/periodicResult.vtk", mySolver->getMesh(), mySolver->getResults(), Vector3(0.0, 0.0551168, 0.0));


    int g = 5;

    return 0;
}