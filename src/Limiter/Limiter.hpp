#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../Mat.hpp"

class Limiter
{
    public:

        Limiter() {}

        virtual ~Limiter() {}

        virtual Field<Vars<5>> calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const;

    protected:
        virtual double limiterFunction(double y) const;

};

#endif // LIMITER_HPP