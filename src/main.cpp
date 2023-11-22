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

    Mesh myMesh = setter.createMesh();
    std::unique_ptr<FluxSolver> myFluxSolver = setter.createFluxSolver();
    std::unique_ptr<Thermo> myThermoModel = setter.createThermoModel();

    std::unique_ptr<FVMScheme> mySolver = setter.createSolver(std::move(myMesh), std::move(myFluxSolver), std::move(myThermoModel));

    mySolver->setCfl(setter.getCfl());
    mySolver->setMaxIter(setter.getMaxIter());
    mySolver->setTargetError(setter.getTargetError());
    mySolver->setLocalTimeStep(setter.getLocalTimeStepSetting());

    std::vector<std::shared_ptr<BoundaryCondition>> bc = setter.createBoundaryCondition(mySolver->getMesh());
    mySolver->setBoundaryConditions(std::move(bc));

    mySolver->setInitialConditions(setter.getInitialCondition());

    outputCFD::outputVTK("../results/results.0.vtk", mySolver->getMesh(), mySolver->getResults());

    mySolver->solve();

    outputCFD::saveFieldOnBoundary("../results/pressure.txt", "wall", mySolver->getMesh(), mySolver->getResults());

    outputCFD::outputVTKPeriodicBoundary("../results/periodicResult.vtk", mySolver->getMesh(), mySolver->getResults(), Vector3(0.0, 0.0551168, 0.0));


    int g = 5;

    return 0;
}