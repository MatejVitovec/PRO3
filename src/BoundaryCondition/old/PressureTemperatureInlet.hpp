#ifndef PRESSURETEMPERATUREINLET_HPP
#define PRESSURETEMPERATUREINLET_HPP

#include "BoundaryCondition.hpp"

class PressureTemperatureInlet : public BoundaryCondition
{
    public:

        PressureTemperatureInlet(Boundary meshBoundary, double totalPressure_, double totalTemperature_, Vars<3> velocityDirection_) : BoundaryCondition(meshBoundary, PRESSURETEMPERATUREINLET),
                    totalPressure(totalPressure_),
                    totalTemperature(totalTemperature_),
                    velocityDirection(velocityDirection_) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;

        void init(const Thermo * const thermoModel);

    private:
        double totalPressure;
        double totalTemperature;
        Vars<3> velocityDirection;

        double totalDensity;
        double totalEntropy;
        double totalEnthalpy;
};

#endif // PRESSURETEMPERATUREINLET