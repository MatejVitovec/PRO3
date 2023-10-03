#ifndef IDEALGAS_HPP
#define IDEALGAS_HPP

#include "EquationOfState.hpp"

class IdealGas : public EquationOfState
{
    public:

        IdealGas() : EquationOfState(), gamma(1.4), R(287.05) {}
        IdealGas(double gamma_, double R_, double cp_) : EquationOfState(), gamma(gamma_), R(R_) {}

        void setGamma(double gamma_);
        void setR(double R_);
        double getGamma() const;
        double getR() const;

        //overwritten virtual
        double pressure(const Compressible& data) const;
        double internalEnergy(const Compressible& data) const;
        double soundSpeed(const Compressible& data) const;

        double density(double totalDensity, double machNumeber2) const;
        double soundSpeed(double pressure, double density) const;
        double machNumber2(double totalPressure, double pressure) const;

        Compressible primitiveToConservative(const Vars<5>& primitive) const;

        Compressible nonLinearSubsonicInletBoundaryState(const Compressible& inDomainState, const Vars<3>& normalVector, double totalPressure, double totaltemperature, Vars<3> velocityDirection) const;

    private:
        double gamma;
        double R;

        double cp() const;

};

#endif // IDEALGAS_HPP