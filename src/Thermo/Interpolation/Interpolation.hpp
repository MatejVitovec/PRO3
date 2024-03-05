#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

class Interpolation
{
    public:
        Interpolation() {}

        virtual ~Interpolation() {}

        virtual double calc(double xx, double yy) const = 0;

        virtual double calcInverseX(double zz, double yy, double guessXX) const = 0;
        virtual double calcInverseY(double xx, double zz, double guessYY) const = 0;
};

#endif // INTERPOLATION_HPP