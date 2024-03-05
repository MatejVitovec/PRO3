#ifndef CUBICLIMITER_HPP
#define CUBICLIMITER_HPP

#include "Limiter.hpp"

class CubicLimiter : public Limiter
{
    public:

        CubicLimiter() : K(1.5), a((K-2)/std::pow(K, 3)), b((3.0 - 2.0*K)/std::pow(K, 2)) {}
        CubicLimiter(double K_) : K(K_), a((K-2)/std::pow(K, 3)), b((3.0 - 2.0*K)/std::pow(K, 2)) {}
        virtual ~CubicLimiter() {}

    private:
        const double K;
        const double a;
        const double b;
        double limiterFunction(double y) const;

};

#endif // CUBICLIMITER_HPP