#ifndef BARTHJESPERSEN_HPP
#define BARTHJESPERSEN_HPP

#include "Limiter.hpp"

class BarthJespersen : public Limiter
{
    public:

        BarthJespersen() {}

        virtual ~BarthJespersen() {}

    private:
        double limiterFunction(double y) const;

};

#endif // BARTHJESPERSEN_HPP