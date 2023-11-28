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

#include "Thermo/IdealGas.hpp"
#include "Thermo/Iapws95Thermo.hpp"
#include "Thermo/Iapws95SpecialGas.hpp"

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

        Compressible getInitialCondition();

        std::unique_ptr<FVMScheme> createSolver();
        std::unique_ptr<FVMScheme> createSolver(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_);

        Mesh createMesh();
        std::unique_ptr<FluxSolver> createFluxSolver();
        std::unique_ptr<Thermo> createThermoModel();
        std::vector<std::shared_ptr<BoundaryCondition>> createBoundaryCondition(const Mesh& mesh);

    private:
        std::vector<std::string> data_;

        void errorMessage(std::string str);

        std::string toLower(std::string str);
        std::string removeWhiteSpaces(std::string str);
        bool stringToBool(std::string str);
        std::vector<std::string> stringArrayToVectorOfStrings(std::string str);
        std::vector<double> vectorStringToDouble(std::vector<std::string> in);
        
        std::string findParameterByKey(std::string key, std::vector<std::string> data);
        std::vector<std::string> findParametersByKey(std::string key, std::vector<std::string> data);

        std::vector<std::string> findObjectNamesInGroup(std::vector<std::string> data);
};



#endif // CASESETTER_HPP