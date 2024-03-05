#ifndef BILINEARINTERPOLATION_HPP
#define BILINEARINTERPOLATION_HPP


#include "Interpolation.hpp"


class BiLinearInterpolation : public Interpolation
{
    public:
    
        BiLinearInterpolation() : n(0), m(0), x(0), y(0), coeffs(0), Interpolation() {}

        BiLinearInterpolation(std::vector<double> xx,
                              std::vector<double> yy,
                              std::function<double(double, double)> f);

        ~BiLinearInterpolation() {}

        double calc(double xx, double yy) const;

    private:

        int n;
        int m;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> dx;
        std::vector<double> dy;

        void calcCoeffs(std::function<double(double, double)> f);

        std::vector<std::array<double, 4>> coeffs;
};

#endif // BILINEARINTERPOLATION_HPP