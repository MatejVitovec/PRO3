#ifndef BILINEARINTERPOLATION_HPP
#define BILINEARINTERPOLATION_HPP


#include "Interpolation.hpp"
//#include "InterpolationCoeffs.hpp"


class BiLinearInterpolation : public Interpolation
{
    public:
    
        BiLinearInterpolation(std::vector<int> xSizes_,
                                             std::vector<int> ySizes_,
                                             std::vector<double> xBoundary_,
                                             std::vector<double> yBoundary_,
                                             TransformationType xTransformation_,
                                             TransformationType yTransformation_) : Interpolation(xSizes_, ySizes_, xBoundary_, yBoundary_, xTransformation_, yTransformation_),
                                                                                    coeffs(std::vector<double>((sizeX*sizeY)*4)) {}

        ~BiLinearInterpolation() {}

        void calcCoeffs(std::function<double(double, double)> f);

        double calc(double x, double y) const;
        double calcInverseX(double y, double z, double xGuess) const;
        double calcInverseY(double x, double z, double yGuess) const;

    private:

        std::vector<double> coeffs;

        std::array<double, 4> getNodeCoeffs(int i, int j) const;

        std::tuple<double, double, double, double> getNodeXYdXdY(int i, int j) const;

};

#endif // BILINEARINTERPOLATION_HPP