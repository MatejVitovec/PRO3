#ifndef PERIODICITY_HPP
#define PERIODICITY_HPP

#include "BoundaryCondition.hpp"

class Periodicity : public BoundaryCondition
{
    public:

        Periodicity() {}
        Periodicity(Boundary meshBoundary, Vector3 faceMidpointShift_, std::string associatedBoundaryName_) : BoundaryCondition(meshBoundary), faceMidpointShift(faceMidpointShift_), associatedBoundaryName(associatedBoundaryName_) {}

        void init(const Mesh& mesh);

        Compressible calculateState(const Compressible& wl, const Face& f, const Thermo * const thermoModel) const;
        void apply(const std::vector<int>& ownerIndexList,const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const;
    
    private:
        std::vector<int> periodicityFacesIndex;
        std::vector<int> periodicityFacesOwnersIndexes;
        std::string associatedBoundaryName; //mozna nepotrebuju - potom zjistit, pripadne odstranit
        Vector3 faceMidpointShift;
};

#endif // PERIODICITY