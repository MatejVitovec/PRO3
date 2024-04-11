#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../Field.hpp"
#include "../Mesh/Mesh.hpp"

#include <chrono>

class Thermo
{
    public:
    
        Thermo() {}

        virtual ~Thermo() {}

        Field<Compressible> updateField(Field<Compressible> w) const;
        Field<Compressible> updateInetrnalFieldFaces(Field<Compressible> w, const Mesh& mesh) const;

        virtual Vars<3> updateThermo(const Compressible& data) const = 0;

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible stagnationState(double TTot, double pTot) const = 0;
        //virtual Compressible isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const = 0;
        //virtual Compressible isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const = 0;

        virtual std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const = 0;

        virtual Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const = 0;
        
};

#endif // THERMO_HPP