#ifndef VARS_HPP
#define VARS_HPP

#include <vector>
#include <memory>
#include <cmath>

#include "Mesh/Vector3.hpp"

template <int N>
class Vars
{
    public:
        Vars() : data() {}
        Vars(const std::array<double, N>& in) : data(in) {}
        Vars(double in) : data() { for (size_t i = 0; i < N; i++) data[i] = in; }

        virtual ~Vars() {}

        inline const double& operator[](int i) const
        {
            return data[i];
        }

        inline double& operator[] (int i)
        {
            return data[i];
        }

        void operator+=(const Vars<N>& v);
        void operator-=(const Vars<N>& v);

        Vars<N> operator-();
        
    protected:
        std::array<double, N> data;
};

template <int N>
void Vars<N>::operator+= (const Vars<N>& v)
{
    for (int i = 0; i < N; i++)
    {
        data[i] += v[i];
    }
}

template <int N>
void Vars<N>::operator-= (const Vars<N>& v)
{
    for (int i = 0; i < N; i++)
    {
        data[i] -= v[i];
    }
}

template <int N>
Vars<N> Vars<N>::operator-()
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = -data[i];
    }
    
    return out;
}

//////////////Non member operators///////////////////

template <int N>
bool operator== (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        if(u[i] != v[i])
            return false;
    }
    return true;
}

template <int N>
Vars<N> operator+ (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i] + v[i];
    }
    return out;
}

// u - v
template <int N>
Vars<N> operator- (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i] - v[i];
    }
    return out;
}

// w * u
template <int N>
Vars<N> operator* (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*v[i];
    }
    return out;
}

// w / u (po slozkach)
template <int N>
Vars<N> operator/ (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]/v[i];
    }
    return out;
}

// a * u
template <int N>
Vars<N> operator* (double a, const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*a;
    }
    return out;
}

// u * a
template <int N>
Vars<N> operator* (const Vars<N>& u, double a)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*a;
    }
    return out;
}

// u / a
template <int N>
Vars<N> operator/ (const Vars<N>& u, double a)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]/a;
    }
    return out;
}


//////////////Non member function///////////////////

template <int N>
Vars<N> sqrt(const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = sqrt(u[i]);
    }
    return out;
}

template <int N>
Vars<N> abs(const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = fabs(u[i]);
    }
    return out;
}

template <int N>
double dot(const Vars<N>& u, const Vars<N>& v)
{
    double out = 0;
    for (int i = 0; i < N; i++)
    {
        out += u[i]*v[i];
    }
    return out;
}

template <int N>
double sum(const Vars<N>& u)
{
    double out;
    for (int i = 0; i < N; i++)
    {
        out += u[i];
    }
    return out;
}

template <int N>
Vars<N> max(const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = std::max(u[i], v[i]);
    }
    return out;
}

template <int N>
Vars<N> min(const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = std::min(u[i], v[i]);
    }
    return out;
}

template <int N>
Vars<N> zeroSmallNumbers(const Vars<N>& u)
{
    constexpr double SMALLNUM = 1e-12;

    Vars<N> out = u;
    for (int i = 0; i < N; i++)
    {
        if(std::fabs(out[i]) < SMALLNUM)
        {
            out[i] = 0.0;
        }
    }
    return out;
}

Vars<3> vector3toVars(const Vector3& vec);


#endif // VARS_HPP