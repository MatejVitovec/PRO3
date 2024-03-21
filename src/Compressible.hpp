#ifndef COMPRESSIBLE_HPP
#define COMPRESSIBLE_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"

class Compressible : public Vars<5>
{
    public:
        enum {RHO, RHO_U, RHO_V, RHO_W, RHO_E};

        Compressible() : Vars<5>() {}
        Compressible(const std::array<double, 5>& in) : Vars<5>(in) {}
        Compressible(const std::array<double, 5>& in, const std::array<double, 3>& inTermo) : Vars<5>(in), thermoVar(inTermo) {}

        void setThermoVar(Vars<3> thermoProp); //TODO prejmenovat na setThernm

        virtual ~Compressible() {}

        void operator+=(const Compressible& v);
        void operator-=(const Compressible& v);
        void operator+=(const Vars<5>& v);
        void operator-=(const Vars<5>& v);

        double density() const;
        Vars<3> velocity() const;
        double absVelocity() const;
        double absVelocity2() const;
        double normalVelocity(const Vars<3>& normalVector) const;
        double velocityU() const;
        double velocityV() const;
        double velocityW() const;
        double totalEnergy() const;
        double internalEnergy() const;

        double temperature() const;
        double pressure() const;        
        double soundSpeed() const;
        double machNumber() const;
        Vars<3> thermo() const;

        Vars<5> flux(const Vars<3>& normalVector) const;
        Vars<5> primitive() const;

    private:
        enum {T, P, A};
        Vars<3> thermoVar;
};



//////////////Non member operators///////////////////

// u == v
bool operator== (const Compressible& u, const Compressible& v);

// u + v
Compressible operator+ (const Compressible& u, const Compressible& v);

// u - v
Compressible operator- (const Compressible& u, const Compressible& v);

// u * v
Compressible operator* (const Compressible& u, const Compressible& v);

// u / v
Compressible operator/ (const Compressible& u, const Compressible& v);

// a * u
Compressible operator* (double a, const Compressible& u);

// u * a
Compressible operator* (const Compressible& u, double a);

// u / a
Compressible operator/ (const Compressible& u, double a);

// Compressible, Vars

// u + v
Compressible operator+ (const Compressible& u, const Vars<5>& v);

// u - v
Compressible operator- (const Compressible& u, const Vars<5>& v);

// w * u
Compressible operator* (const Compressible& u, const Vars<5>& v);



//////////////Non member function///////////////////

Compressible sqrt(const Compressible& u);

#endif // COMPRESSIBLE_HPP