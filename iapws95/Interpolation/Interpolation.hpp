#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <functional>
#include <utility>


class Interpolation
{
    public:

        enum TransformationType{NONE, LOG, LOG10, LOGINV};

        Interpolation(std::vector<int> xSizes_,
                      std::vector<int> ySizes_,
                      std::vector<double> xBoundary_,
                      std::vector<double> yBoundary_,
                      TransformationType xTransformation_,
                      TransformationType yTransformation_);

        virtual ~Interpolation() {}

        double operator()(double x, double y) const
        {
            return calc(x, y);
        }

        virtual void calcCoeffs(std::function<double(double, double)> f) = 0;
        virtual double calc(double x, double y) const = 0;
        virtual double calcInverseX(double y, double z, double xGuess) const = 0;
        virtual double calcInverseY(double x, double z, double yGuess) const = 0;

    protected:

        TransformationType xTransformation;
        TransformationType yTransformation;

        std::vector<int> xSizes;
        std::vector<int> ySizes;
        std::vector<double> xBoundary;
        std::vector<double> yBoundary;

        int sizeX;
        int sizeY;

        std::vector<double> dx;
        std::vector<double> dy;


        double xMax() const;
        double yMax() const;

        double transformX(double x) const;
        double transformY(double y) const;
        double backTransformX(double x) const;
        double backTransformY(double y) const;

        bool checkInterval(double x, double y) const;

        std::pair<int, int> findPosition(double x, double y) const;
};

#endif // INTERPOLATION_HPP