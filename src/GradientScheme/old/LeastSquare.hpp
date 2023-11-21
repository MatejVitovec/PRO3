#ifndef LEASTSQUARE_HPP
#define LEASTSQUARE_HPP

#include "GradientScheme.hpp"

class LeastSquare : public GradientScheme
{
    public:

        LeastSquare() {}

        virtual ~LeastSquare() {}

        void init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

        Field<std::array<Vars<3>, 5>> calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const;

    private:
        Field<std::array<Vars<3>, 3>> MInv;

        void calculateInverseM(Field<std::array<Vars<3>, 3>> M);        
        double det3by3(std::array<Vars<3>, 3> M) const;
        std::array<Vars<3>, 3> adjoint3by3(std::array<Vars<3>, 3> M) const;

        

};

#endif // LEASTSQUARE_HPP