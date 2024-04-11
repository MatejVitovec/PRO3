#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "../Mesh/Mesh.hpp"

#include "../Mesh/Boundary.hpp"
#include "../Mesh/Face.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../Mat.hpp"
#include "../Thermo/Thermo.hpp"

class BoundaryCondition
{
    public:
        enum BoundaryConditionType{ISENTROPICINLET, PRESSUREOUTLET, FREEBOUNDARY, WALL, PERIODICITY};

        BoundaryCondition(BoundaryConditionType type_) : type(type_) {}
        BoundaryCondition(Boundary meshBoundary, BoundaryConditionType type_) : boundary(meshBoundary), type(type_) {}

        BoundaryConditionType getType() const;
        Boundary getBoundary() const;

        //TODO predelat na mesh, wr na konec argumentu funkce
        virtual void apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const;

        virtual void correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;

        virtual Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const = 0;
        
    protected:
        Boundary boundary;
        BoundaryConditionType type;
};

#endif // BOUNDARYCONDITION_HPP