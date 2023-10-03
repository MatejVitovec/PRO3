#ifndef PRESSUREOUTLET_HPP
#define PRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class PressureOutlet : public BoundaryCondition
{
    public:

        PressureOutlet() {}
        PressureOutlet(Boundary meshBoundary, double pressure_) : BoundaryCondition(meshBoundary), pressure(pressure_) {}

        void setPressure(double pressure_);
        double getPressure() const;

        Compressible calculateState(const Compressible& wl, const Face& f) const;
        

    private:
        double pressure;
};

#endif // PRESSUREOUTLET