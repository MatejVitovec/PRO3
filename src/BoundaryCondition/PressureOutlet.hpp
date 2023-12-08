#ifndef PRESSUREOUTLET_HPP
#define PRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class PressureOutlet : public BoundaryCondition
{
    public:

        PressureOutlet(Boundary meshBoundary, double pressure_) : BoundaryCondition(meshBoundary, PRESSUREOUTLET), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;
        

    private:
        double pressure;
};

#endif // PRESSUREOUTLET