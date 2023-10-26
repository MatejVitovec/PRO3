#ifndef LEASTSQUARE_HPP
#define LEASTSQUARE_HPP

#include "GradientScheme.hpp"

class LeastSquare : public GradientScheme
{
    public:

        LeastSquare() {}

        virtual ~LeastSquare() {}

        void init(const Mesh& mesh);

        Field<std::array<Vars<5>, 3>> calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const;

    private:
        Field<Vars<3>> xRowOfMInv;
        Field<Vars<3>> yRowOfMInv;
        Field<Vars<3>> zRowOfMInv;

        Field<Vars<3>> delta;

        void calculatesDeltas(const Mesh& mesh);

        void calculateInverseM(Field<Vars<3>> xRowOfM, Field<Vars<3>> yRowOfM, Field<Vars<3>> zRowOfM);        
        double det3by3(Vars<3> xRow, Vars<3> yRow, Vars<3> zRow) const;
        std::array<double, 9> adjoint3by3(Vars<3> xRow, Vars<3> yRow, Vars<3> zRow) const;

};

#endif // LEASTSQUARE_HPP