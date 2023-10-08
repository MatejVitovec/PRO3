#include "BoundaryCondition.hpp"




void BoundaryCondition::apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const
{
    for (auto & faceIndex : boundary.facesIndex)
    {
        wr[faceIndex] = calculateState(w[ownerIndexList[faceIndex]], faces[faceIndex], thermoModel);
        //wr[faceIndex].setThermoVar(thermoModel->updateThermo(wr[faceIndex])); //zatim docasne (presun do calculateState-konkretni BC - optimalizace)
    }    
}