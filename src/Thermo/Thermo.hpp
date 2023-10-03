#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../Field.hpp"


class Thermo
{
    public:
    
        Thermo() {}

        virtual Vars<3> calculateThermo(const Compressible& data) const = 0;
        Field<Compressible> updateField(const Field<Compressible>& w) const;

};

#endif // THERMO_HPP