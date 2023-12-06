#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../Field.hpp"
#include "../Mesh/Mesh.hpp"


class Thermo
{
    public:
    
        Thermo() {}

        virtual ~Thermo() {}

        Field<Compressible> updateField(Field<Compressible> wn, const Field<Compressible>& w) const;
        Field<Compressible> updateInetrnalFieldFaces(Field<Compressible> wn, const Field<Compressible>& w, const Mesh& mesh) const;

        virtual Vars<3> updateThermo(const Compressible& data, const Compressible& dataOld) const = 0;
        virtual Vars<3> updateThermo(const Compressible& data) const = 0;

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible stagnationState(double TTot, double pTot) const = 0;
        virtual Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const = 0;
        virtual Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const = 0;
        
};

#endif // THERMO_HPP