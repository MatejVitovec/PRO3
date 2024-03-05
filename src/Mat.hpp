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
        Mat() : data_() {}

        virtual ~Mat() {}

        double* operator[](int row)
        {
            return data_.data() + row*N;
        }

        const double* operator[](int row) const
        {
            return data_.data() + row*N;
        }

        /*double* data()
        {
            return data_.data();
        }*/

        void operator+=(const Mat<M, N>& v);
        void operator-=(const Mat<M, N>& v);
        
    protected:
        std::array<double, M*N> data_;
};

template <int M, int N>
void Mat<M, N>::operator+= (const Mat<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data_[i] += v.data_[i];
    }
}

template <int M, int N>
void Mat<M, N>::operator-= (const Mat<M, N>& v)
{
    for (int i = 0; i < M*N; i++)
    {
        data_[i] -= v.data_[i];
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
                out[i][k] += u[i][j]*v[j][k];
            }
        }
    }

    return out;
}

template <int M, int N>
Vars<M> dot(const Mat<M, N>& u, const Vars<N>& v)
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

template <int M, int N>
Mat<M, N> zeroSmallNumbers(const Mat<M, N>& u)
{
    constexpr double SMALLNUM = 1e-15;

    Mat<M, N> out = u;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(std::fabs(out[i][j]) < SMALLNUM)
            {
                out[i][j] = 0.0;
            }
        }
    }
    return out;
}

///////// 3X3

Mat<3,3> adj(const Mat<3,3>& u);

double det(const Mat<3,3>& u);

Mat<3, 3> inv(const Mat<3, 3>& u);

#endif // MAT_HPPX