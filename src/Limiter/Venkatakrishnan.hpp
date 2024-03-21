#ifndef VENKATAKRISHNAN_HPP
#define VENKATAKRISHNAN_HPP

#include "Limiter.hpp"

class Venkatakrishnan : public Limiter
{
    public:

        Venkatakrishnan() {}

        virtual ~Venkatakrishnan() {}

    private:
        double limiterFunction(double y) const;

};

#endif // VENKATAKRISHNAN_HPP