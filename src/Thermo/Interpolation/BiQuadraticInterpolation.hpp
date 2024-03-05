#ifndef BIQUADRATICINTERPOLATION_HPP
#define BIQUADRATICINTERPOLATION_HPP

#include <vector>
#include <functional>

#include "Interpolation.hpp"
#include "Matrix.hpp"

class BiQuadraticInterpolation : public Interpolation
{
    public:

        BiQuadraticInterpolation() : n(0), m(0), x(0), y(0), coeffs(0), Interpolation() {}

        BiQuadraticInterpolation(std::vector<double> xx,
                                 std::vector<double> yy,
                                 std::function<double(double, double)> f);

        BiQuadraticInterpolation(std::vector<double> xx,
                                 std::vector<double> yy,
                                 std::function<double(double, double)> f,
                                 std::function<double(double, double)> fx,
                                 std::function<double(double, double)> fy,
                                 std::function<double(double, double)> fxy);
        
        ~BiQuadraticInterpolation() {}

        double calc(double xx, double yy) const;

    private:
        class Mat3x3
        {
            public:
                Mat3x3() : data_() {}
                Mat3x3(std::array<double, 9> data) : data_(data) {}
                ~Mat3x3() {}

                double* operator[](int i)
                {
                    return data_.data() + i*3;
                }

                const double* operator[](int i) const
                {
                    return data_.data() + i*3;
                }

                Mat3x3 operator*(Mat3x3 v)
                {
                    Mat3x3 out = Mat3x3();
                    for (int i = 0; i < 3; i++) {
                        for (int k = 0; k < 3; k++) {
                            for (int j = 0; j < 3; j++) { out[i][k] += (data_.data() + i*3)[j]*v[j][k]; }
                        }
                    }
                    return out;
                }

                Mat3x3 transpose()
                {
                    Mat3x3 out = Mat3x3();
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++) { out[i][j] = (data_.data() + j*3)[i]; }
                    }
                    return out;
                }

            private:
                std::array<double, 9> data_;
        };

        int n;
        int m;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> t;
        std::vector<double> dx;
        std::vector<double> dy;

        std::vector<Mat3x3> coeffs;

        void calcCoeffs(std::function<double(double, double)> f,
                        std::function<double(double, double)> fx,
                        std::function<double(double, double)> fy,
                        std::function<double(double, double)> fxy);

        void calcCoeffs(std::function<double(double, double)> f);

        std::vector<double> solveTridiagonal(const std::vector<double>& L, const std::vector<double>& D, const std::vector<double>& U, const std::vector<double>& b) const;
        double calcDerivativeX(std::function<double(double, double)> f, double x, double y, double step) const;
        double calcDerivativeY(std::function<double(double, double)> f, double x, double y, double step) const;
        double calcDerivativeXY(std::function<double(double, double)> f, double x, double y, double stepX, double stepY) const;

};


#endif // BiQuadraticInterpolation_HPP