#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>
#include <sstream>

struct Vector3
{
    Vector3() : x(0.0), y(0.0), z(0.0) {};
    Vector3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {};

    double x;
    double y;
    double z;
};

// u + v
inline Vector3 operator+ (const Vector3& u, const Vector3& v)
{
    return Vector3(u.x + v.x, u.y + v.y, u.z + v.z);
}

// u - v
inline Vector3 operator- (const Vector3& u, const Vector3& v)
{
    return Vector3(u.x - v.x, u.y - v.y, u.z - v.z);
}

// a * u
inline Vector3 operator* (double a, const Vector3& u)
{
    return Vector3(a*u.x, a*u.y, a*u.z);
}

// w * u
inline Vector3 operator* (const Vector3& u, const Vector3& v)
{
    return Vector3(u.x*v.x, u.y*v.y, u.z*v.z);
}

// u / a
inline Vector3 operator/ (const Vector3& u, double a)
{
    return Vector3(u.x/a, u.y/a, u.z/a);
}

//comparison operator
/*bool operator== (const Vector3& u, const Vector3& v)
{
    return (u.x == v.x && u.y == v.y && u.z == v.z);
}*/

//dot product
inline double dot(const Vector3& u, const Vector3& v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

//cross product
inline Vector3 cross(const Vector3& u, const Vector3& v)
{
    return Vector3(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}

//euclidian norm
inline double norm2(const Vector3& u)
{
    return sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
}

//unit vector
inline Vector3 unit(const Vector3& u)
{
    double norm = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    return Vector3(u.x/norm, u.y/norm, u.z/norm);
}

//absolut value of x, y, z
inline Vector3 abs(const Vector3& u)
{
    return Vector3(fabs(u.x), fabs(u.y), fabs(u.z));
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& u)
{
    os << u.x << " " << u.y << " " << u.z;
    return os;
};

#endif // VECTOR3_HPP