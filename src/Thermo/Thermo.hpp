#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../Field.hpp"


class Thermo
{
    public:
    
        Thermo() {}

        Field<Compressible> updateField(Field<Compressible> w) const;

        virtual Vars<3> updateThermo(const Compressible& data) const = 0;

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible isentropicInlet(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const = 0;
        
};

#endif // THERMO_HPP