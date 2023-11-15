#ifndef MAT_HPP
#define MAT_HPP

#include <vector>
#include <memory>
#include <cmath>

#include "Vars.hpp"

template <int M, int N>
class Mat
{
    public:
        Mat() : data() {}

        virtual ~Mat() {}

        double* operator[](int row)
        {
            return data.data() + row*N;
        }

        const double* operator[](int row) const
        {
            return data.data() + row*N;
        }

        double* data()
        {
            return data.data();
        }

        void operator+=(const Mat<M, N>& v);
        void operator-=(const Mat<M, N>& v);
        
    protected:
        std::array<double, N*M> data;
};

template <int M, int N>
void Mat<M, N>::operator+= (const Mat<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data[i] += v.data()[i];
    }
}

template <int M, int N>
void Mat<M, N>::operator-= (const Mat<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data[i] -= v.data()[i];
    }
}
//////////////Non member operators///////////////////

template <int M, int N>
Mat<M, N> operator+ (const Mat<M, N>& u, const Mat<M, N>& v)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j] + v[i][j];
        }
    }
    return out;
}

// u - v
template <int M, int N>
Mat<M, N> operator- (const Mat<M, N>& u, const Mat<M, N>& v)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j] - v[i][j];
        }
    }
    return out;
}

// u * v (po slozkach)
template <int M, int N>
Mat<M, N> operator* (const Mat<M, N>& u, const Mat<M, N>& v)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j]*v[i][j];
        }
    }
    return out;
}

// u * v (po slozkach)
template <int M, int N>
Mat<M, N> operator/ (const Mat<M, N>& u, const Mat<M, N>& v)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j]/v[i][j];
        }
    }
    return out;
}

// a * u
template <int M, int N>
Mat<M, N> operator* (double a, const Mat<M, N>& u)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j]*a;
        }
    }
    return out;
}

// u * a
template <int M, int N>
Mat<M, N> operator* (const Mat<M, N>& u, double a)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j]*a;
        }
    }
    return out;
}

// u / a
template <int M, int N>
Mat<M, N> operator/ (const Mat<M, N>& u, double a)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i][j]/a;
        }
    }
    return out;
}


//////////////Non member function///////////////////

template <int M, int N, int K>
Mat<M, K> dot(const Mat<M, N>& u, const Mat<N, K>& v)
{
    Mat<M, K> out;

    for (int i = 0; i < M; i++)
    {
        for (int k = 0; k < K; k++)
        {
            for (int j = 0; j < N; j++)
            {
                out[i][k] += u[i][n]*v[n][k];
            }
        }
    }

    return out;
}

template <int M, int N>
Vars<M> dot(const Mat<M, N>& u, const Vars<M>& v)
{
    Vars<M> out;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i] += u[i][j]*v[j];
        }
    }

    return out;
}

template <int M, int N>
Mat<M, N> outerProd(const Vars<M>& u, const Vars<N>& v)
{
    Mat<M, N> out;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[i]*v[j];
        }
    }

    return out;
}

template <int M, int N>
Mat<M, N> transpose(const Mat<N, M>& u)
{
    Mat<M, N> out;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = u[j][i];
        }
    }
    return out;
}


///////// 3X3

Mat<3, 3> inv(const Mat<3, 3>& u)
{
    return adj(u)/det(u);
}

Mat<3,3> adj(const Mat<3,3>& u)
{
    Mat<3,3> out;

    out[0][0] = u[1][1]*u[2][2] - u[1][2]*u[2][1];        
    out[0][1] = -(u[1][0]*u[2][2] - u[1][2]*u[2][0]);
    out[0][2] = u[1][0]*u[2][1] - u[1][1]*u[2][0];
    out[1][0] = -(u[0][1]*u[2][2] - u[0][2]*u[2][1]);
    out[1][1] = u[0][0]*u[2][2] - u[0][2]*u[2][0];
    out[1][2] = -(u[0][0]*u[2][1] - u[0][1]*u[2][0]);
    out[2][0] = u[0][1]*u[1][2] - u[0][2]*u[1][1];
    out[2][1] = -(u[0][0]*u[1][2] - u[0][2]*u[1][0]);
    out[2][2] = u[0][0]*u[1][1] - u[0][1]*u[1][0];

    return out;
}

double det(const Mat<3,3>& u)
{
    return u[0][0]*u[1][1]*u[2][2] + u[0][1]*u[1][2]*u[2][0] + u[0][2]*u[1][0]*u[2][1] - u[0][0]*u[1][2]*u[2][1] - u[0][1]*u[1][0]*u[2][2] - u[0][2]*u[1][1]*u[2][0];
}

#endif // MAT_HPPX