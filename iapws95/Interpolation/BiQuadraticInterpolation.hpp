#ifndef BIQUADRATICINTERPOLATION_HPP
#define BIQUADRATICINTERPOLATION_HPP


#include "Interpolation.hpp"
//#include "InterpolationCoeffs.hpp"


class BiQuadraticInterpolation : public Interpolation
{
    public:
    
        BiQuadraticInterpolation(std::vector<int> xSizes_,
                                             std::vector<int> ySizes_,
                                             std::vector<double> xBoundary_,
                                             std::vector<double> yBoundary_,
                                             TransformationType xTransformation_,
                                             TransformationType yTransformation_) : Interpolation(xSizes_, ySizes_, xBoundary_, yBoundary_, xTransformation_, yTransformation_),
                                                                                    coeffs(std::vector<BiQuadraticInterpolationCoeffs>(sizeX*sizeY)) {}

        ~BiQuadraticInterpolation() {}

        void calcCoeffs(std::function<double(double, double)> f);
        
        void calcCoeffs(std::function<double(double, double)> f,
                        std::function<double(double, double)> fx,
                        std::function<double(double, double)> fy,
                        std::function<double(double, double)> fxy);

        double calc(double x, double y) const;
        double calcInverseX(double y, double z, double xGuess) const;
        double calcInverseY(double x, double z, double yGuess) const;

    private:

        class BiQuadraticInterpolationCoeffs
        {
            public:
                double* operator[](int row)
                {
                    return data_.data() + row*3;
                }

                const double* operator[](int row) const
                {
                    return data_.data() + row*3;
                }
        
            protected:
                std::array<double, 9> data_;
        };
        

        std::vector<BiQuadraticInterpolationCoeffs> coeffs;

        BiQuadraticInterpolationCoeffs getNodeCoeffs(int i, int j) const;
        std::pair<double, double> getNodeXY(int i, int j) const;

        std::vector<double> solveTridiagonalSystem(const std::vector<double>& L, const std::vector<double>& D, const std::vector<double>& U, const std::vector<double>& b);
};

#endif // BIQUADRATICINTERPOLATION_HPP