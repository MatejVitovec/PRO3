#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "../Mesh/Mesh.hpp"

#include "../Mesh/Boundary.hpp"
#include "../Mesh/Face.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../Thermo/Thermo.hpp"

class BoundaryCondition
{
    public:
        enum BoundaryConditionType{PRESSURETEMPERATUREINLET, PRESSUREDENSITYINLET, PRESSUREOUTLET, FREEBOUNDARY, WALL, PERIODICITY};

        BoundaryCondition() {}
        BoundaryCondition(Boundary meshBoundary) : boundary(meshBoundary) {}

        virtual void apply(const std::vector<int>& ownerIndexList,const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const;

        virtual Compressible calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const = 0;
        
    protected:
        Boundary boundary;
};

#endif // BOUNDARYCONDITION_HPP