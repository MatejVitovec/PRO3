#ifndef PRESSUREOUTLET_HPP
#define PRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class PressureOutlet : public BoundaryCondition
{
    public:

        PressureOutlet(Boundary meshBoundary, double pressure_) : BoundaryCondition(meshBoundary, PRESSUREOUTLET), pressure(pressure_) {}

        void setPressure(double pressure_);
        double getPressure() const;

        Compressible calculateState(const Compressible& wl, const Compressible& wr, const Face& f, const Thermo * const thermoModel) const;
        

    private:
        double pressure;
};

#endif // PRESSUREOUTLET