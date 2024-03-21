#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP

#include <vector>
#include <memory>

#include "Compressible.hpp"

class Primitive : public Vars<5>
{
    public:
        enum {RHO, U, V, W, P};

        Primitive() : Vars<5>() {}
        Primitive(const std::array<double, 5>& in) : Vars<5>(in) {}
        Primitive(const std::array<double, 5>& in, const std::array<double, 3>& inTermo) : Vars<5>(in), thermoVar(inTermo) {}
        Primitive(Compressible in);

        void setThermo(Vars<3> thermoProp);

        virtual ~Primitive() {}

        void operator+=(const Primitive& v);
        void operator-=(const Primitive& v);
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

        double temperature() const;
        double pressure() const;        
        double soundSpeed() const;
        double machNumber() const;
        double totalEnergy() const;

        Vars<5> conservativeFlux(const Vars<3>& normalVector) const;

    private:
        enum {T, E, A};
        Vars<3> thermoVar;
};



//////////////Non member operators///////////////////

// u == v
bool operator== (const Primitive& u, const Primitive& v);

// u + v
Primitive operator+ (const Primitive& u, const Primitive& v);

// u - v
Primitive operator- (const Primitive& u, const Primitive& v);

// w * u
Primitive operator* (const Primitive& u, const Primitive& v);

// a * u
Primitive operator* (double a, const Primitive& u);

// u * a
Primitive operator* (const Primitive& u, double a);

// u / a
Primitive operator/ (const Primitive& u, double a);

// Primitive, Vars

// u + v
Primitive operator+ (const Primitive& u, const Vars<5>& v);

// u - v
Primitive operator- (const Primitive& u, const Vars<5>& v);

// w * u
Primitive operator* (const Primitive& u, const Vars<5>& v);



//////////////Non member function///////////////////

Primitive sqrt(const Primitive& u);

#endif // PRIMITIVE_HPP