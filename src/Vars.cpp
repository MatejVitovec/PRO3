#include "Vars.hpp"

Vars<3> vector3toVars(const Vector3& vec)
{
    return Vars<3>({vec.x, vec.y, vec.z});
}