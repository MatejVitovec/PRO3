#ifndef ISENTROPICINLET_HPP
#define ISENTROPICINLET_HPP

#include "BoundaryCondition.hpp"

class IsentropicInlet : public BoundaryCondition
{
    public:

        IsentropicInlet(Boundary meshBoundary,
                        double totalDensity_,
                        double totalPressure_,
                        double totalTemperature_,
                        double totalEntropy_,
                        double totalEnthalpy_,
                        Vars<3> velocityDirection_) : BoundaryCondition(meshBoundary, ISENTROPICINLET),
                                                      totalDensity(totalDensity_),
                                                      totalPressure(totalPressure_),
                                                      totalTemperature(totalTemperature_),
                                                      totalEntropy(totalEntropy_),
                                                      totalEnthalpy(totalEnthalpy_),
                                                      velocityDirection(velocityDirection_) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;

    private:
        double totalDensity;
        double totalPressure;
        double totalTemperature;        
        double totalEntropy;
        double totalEnthalpy;
        Vars<3> velocityDirection;
};

#endif // ISENTROPICINLET