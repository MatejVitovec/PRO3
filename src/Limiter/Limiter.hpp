#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../Primitive.hpp"
#include "../Mat.hpp"

class Limiter
{
    public:

        Limiter() {}

        virtual ~Limiter() {}

        virtual Field<Vars<5>> calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const;
        virtual Field<Vars<5>> calculateLimiter(const Field<Primitive>& ul, const Field<Primitive>& ur, const Field<Mat<5,3>>& grad, const Mesh& mesh) const;


    protected:

};

#endif // LIMITER_HPP