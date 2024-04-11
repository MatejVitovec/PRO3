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

        BiLinearInterpolation(std::vector<int> gridSizeX_, std::vector<int> gridSizeY_,
                              std::vector<double> boundaryX_, std::vector<double> boundaryY_,
                              Transformation transformationX_, Transformation transformationY_,
                              std::function<double(double, double)> f);

        ~BiLinearInterpolation() {}

        std::pair<int, int> findPosition(double xx, double yy) const;
        std::pair<int, int> fastFindPosition(double xx, double yy) const;

        double calc(double xx, double yy) const;
        double calcFastFind(double xx, double yy) const;

        double calcFastFindDiffX(double xx, double yy) const;
        double calcFastFindDiffY(double xx, double yy) const;

        double calcInverseX(double zz, double yy, double guessXX) const;
        double calcInverseY(double xx, double zz, double guessYY) const;

    private:

        int n;
        int m;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> dx;
        std::vector<double> dy;

        std::vector<double> dxTransf;
        std::vector<double> dyTransf;

        void calcCoeffs(std::function<double(double, double)> f);

        std::vector<std::array<double, 4>> coeffs;
};

#endif // BILINEARINTERPOLATION_HPP