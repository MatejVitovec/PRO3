#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

class Interpolation
{
    public:
        Interpolation() {}

        virtual ~Interpolation() {}

        virtual double calc(double xx, double yy) const = 0;
};

#endif // INTERPOLATION_HPP