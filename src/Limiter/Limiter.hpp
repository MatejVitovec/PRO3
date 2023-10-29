#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

class Limiter
{
    public:

        Limiter() {}

        virtual Field<double> calculateLimiter(const Mesh& mesh) const;


    protected:

};

#endif // LIMITER_HPP