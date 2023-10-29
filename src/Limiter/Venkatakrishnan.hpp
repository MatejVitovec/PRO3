#ifndef VENKATAKRISHNAN_HPP
#define VENKATAKRISHNAN_HPP

#include "Limiter.hpp"

class Venkatakrishnan : public Limiter
{
    public:

        Venkatakrishnan() {}

        virtual ~Venkatakrishnan() {}

        Field<double> calculateLimiter(const Mesh& mesh) const;


    private:


};

#endif // VENKATAKRISHNAN_HPP