#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>

#include "Mesh/Mesh.hpp"
#include "ExplicitEuler.hpp"
#include "FluxSolver/Hll.hpp"
#include "FluxSolver/Hllc.hpp"
#include "Thermo/IdealGas.hpp"
#include "Thermo/Iapws95.hpp"
#include "outputCFD.hpp"
#include "setCFD.hpp"



int main(int argc, char** argv)
{
    Mesh myMesh = Mesh();

    //auto stop1 = std::chrono::high_resolution_clock::now();

    //myMesh.loadGmsh2("../meshes/GAMM.msh");
    myMesh.loadGmsh2("../meshes/nozzle.msh");

    /*auto stop2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1).count() << " ms\n";*/


    std::unique_ptr<FluxSolver> myFluxSolver = std::make_unique<Hllc>();
    //std::unique_ptr<Thermo> myThermoModel = std::make_unique<IdealGas>(1.4, 461.51805);
    std::unique_ptr<Thermo> myThermoModel = std::make_unique<Iapws95>();

    ExplicitEuler mySolver(std::move(myMesh), std::move(myFluxSolver), std::move(myThermoModel));

    mySolver.setCfl(0.6);
    mySolver.setMaxIter(5000000);
    mySolver.setTargetError(0.0000005);
    mySolver.setLocalTimeStep(false);

    std::vector<std::unique_ptr<BoundaryCondition>> bc = createBoundaryCondition(mySolver.getMesh());

    mySolver.setBoundaryConditions(std::move(bc));

    mySolver.setInitialConditionsPrimitive(Vars<5>({0.5, 0.0, 0.0, 0.0, 80000.0}));

    outputVTK("../results/results.0.vtk", mySolver.getMesh(), mySolver.getResults());

    mySolver.solve();

    int a = 5;

    return 0;
}