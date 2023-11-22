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

        void setTotalPressure(double totalPressure_);
        void setTotalTemperature(double totalTemperature_);
        void setVelocityDirection(Vars<3> velocityDirection_);
        double getTotalPressure() const;
        double getTotalTemperature() const;
        Vars<3> getVelocityDirection() const;

        Compressible calculateState(const Compressible& wl, const Compressible& wrOld, const Face& f, const Thermo * const thermoModel) const;        

    private:
        double totalPressure;
        double totalTemperature;
        Vars<3> velocityDirection;
};

#endif // PRESSURETEMPERATUREINLET