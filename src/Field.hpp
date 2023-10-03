#ifndef FIELD_H
#define FIELD_H

#include <vector>

template <typename T>
class Field
{
    public:
        Field() : data() {}
        Field(int n) : data(n) {}

        int size() const;
        T norm() const;

        T operator[](int i) const
        {
            return data[i];
        }

        T& operator[](int i)
        {
            return data[i];
        }

        template<typename U>
        void operator+=(const Field<U>& v);

        template<typename U>
        void operator-=(const Field<U>& v);

    private:
        std::vector<T> data;
};


template <typename T>
int Field<T>::size() const
{
    return data.size();
}

template <typename T>
T Field<T>::norm() const
{
    T sum;

    for (int i = 0; i < data.size(); i++)
    {
        sum += (data[i]*data[i]);
    }

    return sqrt(sum);
}

template<typename T>
template<typename U>
void Field<T>::operator+=(const Field<U>& v)
{
    for (int i = 0; i < data.size(); i++)
    {
        data[i] += v[i];
    }
}

template<typename T>
template<typename U>
void Field<T>::operator-=(const Field<U>& v)
{
    for (int i = 0; i < data.size(); i++)
    {
        data[i] += v[i];
    }
}


template <typename T, typename U>
Field<T> operator+ (const Field<T>& u, const Field<U>& v)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i] + v[i];
    }

    return out;
}

template <typename T, typename U>
Field<T> operator- (const Field<T>& u, const Field<U>& v)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i] - v[i];
    }

    return out;
}

template <typename T>
Field<T> operator+ (const Field<T>& u, const T& a)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i] - a;
    }

    return out;
}

template <typename T>
Field<T> operator- (const Field<T>& u, const T& a)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i] - a;
    }

    return out;
}

template <typename T, typename D>
Field<T> operator* (const Field<T>& u, const D& a)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i]*a;
    }

    return out;
}

template <typename T, typename D>
Field<T> operator/ (const Field<T>& u, const D& a)
{
    Field<T> out(u.size());

    for (int i = 0; i < u.size(); i++)
    {
        out[i] = u[i]/a;
    }

    return out;
}

#endif // FIELD_H