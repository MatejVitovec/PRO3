#include <cmath>
#include "Primitive.hpp"

Primitive::Primitive(Compressible in)
{
    data[RHO] = in.density();
    data[U] = in.velocityU();
    data[V] = in.velocityV();
    data[W] = in.velocityW();
    data[P] = in.pressure();

    thermoVar[T] = in.temperature();
    thermoVar[A] = in.soundSpeed();
    thermoVar[E] = in.totalEnergy();
}

void Primitive::setThermo(Vars<3> thermoProp)
{
    thermoVar = thermoProp;
}

double Primitive::density() const
{
    return data[RHO];
}

Vars<3> Primitive::velocity() const
{
    return Vars<3>({data[U], data[V], data[W]});
}

double Primitive::absVelocity() const
{
    return sqrt(data[U]*data[U] + data[V]*data[V] + data[W]*data[W]);
}

double Primitive::absVelocity2() const
{
    return data[U]*data[U] + data[V]*data[V] + data[W]*data[W];
}

double Primitive::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double Primitive::velocityU() const
{
    return data[U];
}

double Primitive::velocityV() const
{
    return data[V];
}

double Primitive::velocityW() const
{
    return data[W];
}

double Primitive::temperature() const
{
    return thermoVar[T];
}

double Primitive::pressure() const
{
    return data[P];
}

double Primitive::soundSpeed() const
{
    return thermoVar[A];
}

double Primitive::machNumber() const
{
    return absVelocity()/soundSpeed();
}

double Primitive::totalEnergy() const
{
    return thermoVar[E];
}

Vars<5> Primitive::conservativeFlux(const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return Compressible({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoVar[P] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoVar[P] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoVar[P] * normalVector[2],
                         (data[RHO]*thermoVar[E] + thermoVar[P]) * normalVelocity});
}

void Primitive::operator+=(const Primitive& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void Primitive::operator-=(const Primitive& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}

void Primitive::operator+=(const Vars<5>& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void Primitive::operator-=(const Vars<5>& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}


//////////////Non member operators///////////////////

// u == v
bool operator== (const Primitive& u, const Primitive& v)
{
    if(u[0] == v[0] && u[1] == v[1] && u[2] == v[2] && u[3] == v[3] && u[4] == v[4]) return true;
    return false;
}

// u + v
Primitive operator+ (const Primitive& u, const Primitive& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
Primitive operator- (const Primitive& u, const Primitive& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
Primitive operator* (const Primitive& u, const Primitive& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}

// a * u
Primitive operator* (double a, const Primitive& u)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u * a
Primitive operator* (const Primitive& u, double a)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u / a
Primitive operator/ (const Primitive& u, double a)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] /= a;
    }
    
    return out;
}

// Primitive, Vars<5>

// u + v
Primitive operator+ (const Primitive& u, const Vars<5>& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
Primitive operator- (const Primitive& u, const Vars<5>& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
Primitive operator* (const Primitive& u, const Vars<5>& v)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}


//////////////Non member function///////////////////

Primitive sqrt(const Primitive& u)
{
    Primitive out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] = std::sqrt(u[i]);
    }
    
    return out;
}