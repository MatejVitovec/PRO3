#ifndef VARS2D_HPP
#define VARS2D_HPP

#include <vector>
#include <memory>
#include <cmath>

#include "Mesh/Vector3.hpp"
#include "Vars.hpp"

template <int M, int N>
class Vars2d;

template <int M, int N>
Vars2d<M, N> operator+ (const Vars2d<M, N>& u, const Vars2d<M, N>& v);

template <int M, int N>
Vars2d<M, N> operator- (const Vars2d<M, N>& u, const Vars2d<M, N>& v);

template <int M, int N>
Vars2d<M, N> operator* (const Vars2d<M, N>& u, const Vars2d<M, N>& v);

template <int M, int N>
Vars2d<M, N> operator* (double a, const Vars2d<M, N>& u);

template <int M, int N>
Vars2d<M, N> operator* (const Vars2d<M, N>& u, double a);

template <int M, int N>
Vars2d<M, N> operator/ (const Vars2d<M, N>& u, double a);

template <int M, int N>
Vars2d<M, N> sqrt(const Vars2d<M, N>& u);

template <int M, int N>
Vars2d<M, N> abs(const Vars2d<M, N>& u);

template <int M, int N>
Vars<N> dot(const Vars2d<M, N>& u, const Vars<N>& v);


template <int M, int N>
class Vars2d
{
    public:
        Vars2d() : data() {}

        Vars2d(const std::array<double, N>& in) : data(in) {}

        virtual ~Vars2d() {}

        inline const Vars<N>& operator[](int i) const
        {
            Vars<N> out;
            for (int j = 0; j < N; j++)
            {
                out[j] = data[i*N + j];
            }
            return out;
        }

        inline Vars<N>& operator[] (int i)
        {
            Vars<N> out;
            for (int j = 0; j < N; j++)
            {
                out[j] = data[i*N + j];
            }
            return out;
        }

        void operator+=(const Vars2d<M, N>& v);
        void operator-=(const Vars2d<M, N>& v);

        friend Vars2d<M, N> operator+ <> (const Vars2d<M, N>& u, const Vars2d<M, N>& v);
        friend Vars2d<M, N> operator- <> (const Vars2d<M, N>& u, const Vars2d<M, N>& v);
        friend Vars2d<M, N> operator* <> (const Vars2d<M, N>& u, const Vars2d<M, N>& v);
        friend Vars2d<M, N> operator* <> (double a, const Vars2d<M, N>& u);
        friend Vars2d<M, N> operator* <> (const Vars2d<M, N>& u, double a);
        friend Vars2d<M, N> operator/ <> (const Vars2d<M, N>& u, double a);
        friend Vars2d<M, N> sqrt <> (const Vars2d<M, N>& u);
        friend Vars2d<M, N> abs <> (const Vars2d<M, N>& u);
        friend Vars<N> dot <> (const Vars2d<M, N>& u, const Vars<N>& v);

        
    protected:
        std::array<double, M*N> data;
};

template <int M, int N>
void Vars2d<M, N>::operator+= (const Vars2d<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data[i] += v.data[i];
    }
}

template <int M, int N>
void Vars2d<M, N>::operator-= (const Vars2d<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data[i] -= v.data[i];
    }
}

//////////////Non member operators///////////////////

template <int M, int N>
Vars2d<M, N> operator+ (const Vars2d<M, N>& u, const Vars2d<M, N>& v)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i] + v.data[i];
    }
    return out;
}

// u - v
template <int M, int N>
Vars2d<M, N> operator- (const Vars2d<M, N>& u, const Vars2d<M, N>& v)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i] - v.data[i];
    }
    return out;
}

// w * u
template <int M, int N>
Vars2d<M, N> operator* (const Vars2d<M, N>& u, const Vars2d<M, N>& v)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i] * v.data[i];
    }
    return out;
}

// a * u
template <int M, int N>
Vars2d<M, N> operator* (double a, const Vars2d<M, N>& u)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i]*a;
    }
    return out;
}

// u * a
template <int M, int N>
Vars2d<M, N> operator* (const Vars2d<M, N>& u, double a)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i]*a;
    }
    return out;
}

// u / a
template <int M, int N>
Vars2d<M, N> operator/ (const Vars2d<M, N>& u, double a)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = u.data[i]/a;
    }
    return out;
}

//////////////Non member function///////////////////

template <int M, int N>
Vars2d<M, N> sqrt(const Vars2d<M, N>& u)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = sqrt(u.data[i]);
    }
    return out;
}

template <int M, int N>
Vars2d<M, N> abs(const Vars2d<M, N>& u)
{
    Vars2d<M, N> out;
    for (int i = 0; i < M*N; i++)
    {
        out.data[i] = fabs(u.data[i]);
    }
    return out;
}

template <int M, int N>
Vars<N> dot(const Vars2d<M, N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i] += u[i*N + j]*v[j];
        }
        
    }
    return out;
}


#endif // VARS2D_HPP