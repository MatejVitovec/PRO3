#ifndef CASESETTER_HPP
#define CASESETTER_HPP

#include <vector>
#include <memory>

#include "Mesh/Mesh.hpp"
#include "FluxSolver/FluxSolver.hpp"
#include "Thermo/Thermo.hpp"
#include "Compressible.hpp"
#include "Field.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

#include "ExplicitEuler.hpp"
#include "HeunScheme.hpp"
#include "FluxSolver/Hll.hpp"
#include "FluxSolver/Hllc.hpp"

#include "Thermo/IdealGasThermo.hpp"
#include "Thermo/Iapws95Thermo.hpp"
#include "Thermo/SpecialGasThermo.hpp"
#include "Thermo/Iapws95InterpolationThermo.hpp"
#include "Thermo/Interpolation/BiLinearInterpolation.hpp"
#include "Thermo/Interpolation/BiQuadraticInterpolation.hpp"

#include "GradientScheme/LeastSquare.hpp"
#include "Limiter/BarthJespersen.hpp"
#include "Limiter/Venkatakrishnan.hpp"
#include "Limiter/CubicLimiter.hpp"

/*#include "Thermo/IdealGas.hpp"
#include "Thermo/Iapws95.hpp"
#include "Thermo/Iapws95SpecialGas.hpp"*/


class CaseSetter
{
    public:

        CaseSetter() {}

        void loadSettingFile(std::string fileName);

        bool getLocalTimeStepSetting();
        double getTargetError();
        double getMaxIter();
        double getCfl();
        int getSaveIterInterval();

        Compressible getInitialCondition(const Thermo * const thermoModel);

        std::unique_ptr<FVMScheme> createSolver();
        std::unique_ptr<FVMScheme> createSolver(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_);
        std::unique_ptr<FVMScheme> createAndSetSolver();

        Mesh createMesh();
        std::unique_ptr<FluxSolver> createFluxSolver();
        std::unique_ptr<Thermo> createThermoModel();
        std::unique_ptr<GradientScheme> createReconstructionGradient();
        std::unique_ptr<Limiter> createReconstructionLimiter();

        std::vector<std::shared_ptr<BoundaryCondition>> createBoundaryCondition(const Mesh& mesh, const Thermo * const thermoModel);

    private:
        std::vector<std::string> data_;

        void errorMessage(std::string str);

        std::string toLower(std::string str);
        std::string removeWhiteSpaces(std::string str);
        bool stringToBool(std::string str);
        std::vector<std::string> stringArrayToVectorOfStrings(std::string str);
        std::vector<double> vectorStringToDouble(std::vector<std::string> in);
        
        bool isParameterDefined(std::string key, std::vector<std::string> data);
        std::string findParameterByKey(std::string key, std::vector<std::string> data);
        std::vector<std::string> findParametersByKey(std::string key, std::vector<std::string> data);

        std::vector<std::string> findObjectNamesInGroup(std::vector<std::string> data);
};



#endif // CASESETTER_HPP