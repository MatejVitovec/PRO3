#include <cmath>
#include "Compressible.hpp"

void Compressible::setThermoVar(Vars<3> thermoProp)
{
    thermoVar = thermoProp;
}

double Compressible::density() const
{
    return data[RHO];
}

Vars<3> Compressible::velocity() const
{
    return Vars<3>({data[RHO_U] / data[RHO], data[RHO_V] / data[RHO], data[RHO_W] / data[RHO]});
}

double Compressible::absVelocity() const
{
    return sqrt((data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]));
}

double Compressible::absVelocity2() const
{
    return (data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]);
}

double Compressible::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double Compressible::velocityU() const
{
    return data[RHO_U] / data[RHO];
}

double Compressible::velocityV() const
{
    return data[RHO_V] / data[RHO];
}

double Compressible::velocityW() const
{
    return data[RHO_W] / data[RHO];
}

double Compressible::totalEnergy() const
{
    return data[RHO_E] / data[RHO];
}

double Compressible::temperature() const
{
    return thermoVar[T];
}

double Compressible::pressure() const
{
    return thermoVar[P];
}

double Compressible::internalEnergy() const
{
    return data[RHO_E]/data[RHO] - 0.5*this->absVelocity2();
}

double Compressible::soundSpeed() const
{
    return thermoVar[A];
}

double Compressible::machNumber() const
{
    return absVelocity()/soundSpeed();
}

Vars<3> Compressible::thermo() const
{
    return thermoVar;
}

Vars<5> Compressible::flux(const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return Compressible({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoVar[P] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoVar[P] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoVar[P] * normalVector[2],
                         (data[RHO_E]+ thermoVar[P]) * normalVelocity});  //tady byla entalpie - muze byt blbe
}

Vars<5> Compressible::primitive() const
{
    return Vars<5>({data[RHO],
                    data[RHO_U] / data[RHO],
                    data[RHO_V] / data[RHO],
                    data[RHO_W] / data[RHO],
                    thermoVar[P]});
}

void Compressible::operator+=(const Compressible& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void Compressible::operator-=(const Compressible& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}

void Compressible::operator+=(const Vars<5>& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void Compressible::operator-=(const Vars<5>& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}


//////////////Non member operators///////////////////

// u == v
bool operator== (const Compressible& u, const Compressible& v)
{
    if(u[0] == v[0] && u[1] == v[1] && u[2] == v[2] && u[3] == v[3] && u[4] == v[4]) return true;
    return false;
}

// u + v
Compressible operator+ (const Compressible& u, const Compressible& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
Compressible operator- (const Compressible& u, const Compressible& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
Compressible operator* (const Compressible& u, const Compressible& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}

Compressible operator/ (const Compressible& u, const Compressible& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        if (out[i] != 0)
        {
            out[i] /= v[i];
        }
    }
    
    return out;
}

// a * u
Compressible operator* (double a, const Compressible& u)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u * a
Compressible operator* (const Compressible& u, double a)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u / a
Compressible operator/ (const Compressible& u, double a)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] /= a;
    }
    
    return out;
}

// Compressible, Vars<5>

// u + v
Compressible operator+ (const Compressible& u, const Vars<5>& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
Compressible operator- (const Compressible& u, const Vars<5>& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
Compressible operator* (const Compressible& u, const Vars<5>& v)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}


//////////////Non member function///////////////////

Compressible sqrt(const Compressible& u)
{
    Compressible out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] = std::sqrt(u[i]);
    }
    
    return out;
}