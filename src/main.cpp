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

    /*std::shared_ptr<Iapws95Thermo> thermo = std::make_shared<Iapws95Thermo>();
    for (int i = 0; i < 100; i++)
    {
        std::cout << "M2is: " << thermo->calcM2is(5000000.0, 573.15,  2196600 + i*1) << " p2: " <<  2196600 + i*1 << std::endl;
    }*/

    

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

    return 0;
}