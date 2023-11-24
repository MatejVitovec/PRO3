#include "BoundaryCondition.hpp"


BoundaryCondition::BoundaryConditionType BoundaryCondition::getType() const
{
    return type;
}

Boundary BoundaryCondition::getBoundary() const
{
    return boundary;
}

void BoundaryCondition::apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr,const Field<Compressible>& wrOld, const Thermo * const thermoModel) const
{
    for (auto & faceIndex : boundary.facesIndex)
    {
        wr[faceIndex] = calculateState(w[ownerIndexList[faceIndex]], wrOld[faceIndex], faces[faceIndex], thermoModel);
    }    
}

void BoundaryCondition::correct(const Field<Compressible>& w, Field<Compressible> wl, Field<Compressible>& wr, const Field<Compressible>& wrOld, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{

}