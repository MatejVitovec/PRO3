#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "../Mesh/Boundary.hpp"
#include "../Mesh/Face.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

class BoundaryCondition
{
    public:
        enum BoundaryConditionType{PRESSURETEMPERATUREINLET, PRESSUREDENSITYINLET, PRESSUREOUTLET, FREEBOUNDARY, WALL};

        BoundaryCondition() {}
        BoundaryCondition(Boundary meshBoundary) : boundary(meshBoundary) {}

        void apply(const std::vector<int>& ownerIndexList,const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr) const;

        virtual Compressible calculateState(const Compressible& wl, const Face& f) const = 0;
        
    protected:
        Boundary boundary;
};

#endif // BOUNDARYCONDITION_HPP