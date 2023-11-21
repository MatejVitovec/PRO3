#ifndef VENKATAKRISHNAN_HPP
#define VENKATAKRISHNAN_HPP

#include "Limiter.hpp"

class Venkatakrishnan : public Limiter
{
    public:

        Venkatakrishnan() {}

        virtual ~Venkatakrishnan() {}

        Field<Vars<5>> calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Mesh& mesh) const;


    private:


};

#endif // VENKATAKRISHNAN_HPP