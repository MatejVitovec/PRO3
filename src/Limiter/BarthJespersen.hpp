#ifndef BARTHJESPERSEN_HPP
#define BARTHJESPERSEN_HPP

#include "Limiter.hpp"

class BarthJespersen : public Limiter
{
    public:

        BarthJespersen() {}

        virtual ~BarthJespersen() {}

        Field<Vars<5>> calculateLimiter(const Field<Compressible>& wl, const Field<Compressible>& wr, Field<std::array<Vars<3>, 5>> grad, const Mesh& mesh) const;


    private:


};

#endif // BARTHJESPERSEN_HPP