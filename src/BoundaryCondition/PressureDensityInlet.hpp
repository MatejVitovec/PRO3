#ifndef PRESSUREDENSITYINLET_HPP
#define PRESSUREDENSITYINLET_HPP

#include "BoundaryCondition.hpp"

class PressureDensityInlet : public BoundaryCondition
{
    public:

        PressureDensityInlet() {}
        PressureDensityInlet(Boundary meshBoundary, double totalPressure_, double totaltemperature_, Vars<3> velocityDirection_) : BoundaryCondition(meshBoundary),
                    totalPressure(totalPressure_),
                    totalDensity(totaltemperature_),
                    velocityDirection(velocityDirection_) {}

        void setTotalPressure(double totalPressure_);
        void setTotalDensity(double totaltemperature_);
        void setVelocityDirection(Vars<3> velocityDirection_);
        double getTotalPressure() const;
        double getTotalDensity() const;
        Vars<3> getVelocityDirection() const;

        Compressible calculateState(const Compressible& wl, const Face& f) const;        

    private:
        double totalPressure;
        double totalDensity;
        Vars<3> velocityDirection;
};

#endif // PRESSUREDENSITYINLET