#ifndef FVMSCHEME_HPP
#define FVMSCHEME_HPP

#include <vector>
#include <memory>

#include "Compressible.hpp"
#include "Field.hpp"
#include "Mesh/Mesh.hpp"
#include "FluxSolver/FluxSolver.hpp"
#include "Thermo/Thermo.hpp"
#include "GradientScheme/LeastSquare.hpp"
#include "Limiter/Limiter.hpp"
#include "Limiter/BarthJespersen.hpp"
#include "Limiter/Venkatakrishnan.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

class FVMScheme
{
    public:

        //vznikne obraz meshe - neprekopiruje se
        //FVMScheme() : mesh(Mesh()), w(Field<Compressible>()), wl(Field<Compressible>()), wr(Field<Compressible>()) {}
        FVMScheme(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : mesh(std::move(mesh_)),
                                                                                                            fluxSolver(std::move(fluxSolver_)),
                                                                                                            thermo(std::move(thermo_)),
                                                                                                            w(Field<Compressible>()),
                                                                                                            wl(Field<Compressible>()),
                                                                                                            wr(Field<Compressible>()),
                                                                                                            cfl(0.8), maxIter(10000000),
                                                                                                            targetError(0000005),
                                                                                                            localTimeStep(false),
                                                                                                            reconstruction(false),
                                                                                                            time(0.0),
                                                                                                            gradientScheme(std::make_unique<GradientScheme>()),
                                                                                                            limiter(std::make_unique<Limiter>()) {}


        virtual ~FVMScheme() {}

        void setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_);
        void setReconstructionLimiter(std::unique_ptr<Limiter> limiter_);

        void setCfl(double cfl_);
        void setMaxIter(int maxIter_);
        void setSaveEveryIter(int saveEveryIter_);
        void setTargetError(double targetError_);
        void setLocalTimeStep(bool localTimeStep_);
        void setReconstructionSettings(bool reconstruction_);

        void setThermoModel(std::unique_ptr<Thermo> thermo_) {thermo = std::move(thermo_);}
        
        double getCfl() const;
        int getMaxIter() const;
        double getTargetError() const;
        bool getTimeStepsettings() const;
        bool getReconstructionSettings() const;
        const Mesh& getMesh() const;
        const Thermo* getThermoRef();

        void setInitialConditions(Compressible initialCondition);
        void setInitialConditionsPrimitive(Vars<5> initialCondition);

        void setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions);

        void init();
        virtual void solve() = 0;

        Field<Compressible> getResults() const;
        
        
    protected:
        std::unique_ptr<FluxSolver> fluxSolver;
        std::unique_ptr<Thermo> thermo;

        std::unique_ptr<GradientScheme> gradientScheme;
        std::unique_ptr<Limiter> limiter;

        Mesh mesh;
        std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditionList;

        Field<Compressible> w; //cell size

        Field<Compressible> wl; //faces size
        Field<Compressible> wr;
        Field<Vars<5>> fluxes;

        double cfl;
        int maxIter;
        int saveEveryIter;
        double targetError;
        bool localTimeStep;
        bool reconstruction;

        Field<double> timeSteps;

        double time;

        void updateTimeStep();
        void applyBoundaryConditions();
        void calculateWlWr();
        void calculateFluxes();
        Field<Vars<5>> calculateResidual();
        void reconstruct();
        void boundField();

    private:

};



#endif // FVMSCHEME_HPP