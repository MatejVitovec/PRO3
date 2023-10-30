#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

class Limiter
{
    public:

        Limiter() {}

        virtual ~Limiter() {}

        virtual Field<Vars<5>> calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, Field<std::array<Vars<3>, 5>> grad, const Mesh& mesh) const;


    protected:

};

#endif // LIMITER_HPP