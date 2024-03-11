#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

class Interpolation
{
    public:
        enum Transformation{NONE, LOG, LOG10, LOGINV};
        Interpolation() : gridSizeX(std::vector<int>()), gridSizeY(std::vector<int>()),
                          boundaryX(std::vector<double>()), boundaryY(std::vector<double>()),
                          transformDx(std::vector<double>()), transformDy(std::vector<double>()),
                          transformationX(NONE), transformationY(NONE) {}

        Interpolation(std::vector<int> gridSizeX_, std::vector<int> gridSizeY_,
                      std::vector<double> boundaryX_, std::vector<double> boundaryY_,
                      Transformation transformationX_, Transformation transformationY_) : 
                        gridSizeX(gridSizeX_), gridSizeY(gridSizeY_),
                        boundaryX(boundaryX_), boundaryY(boundaryY_),
                        transformDx(transformationX_), transformDy(transformationY_),
                        transformationX(transformationX_), transformationY(transformationY_) {}

        virtual ~Interpolation() {}

        virtual double calc(double xx, double yy) const = 0;
        virtual double calcFastFind(double xx, double yy) const = 0;

        virtual double calcInverseX(double zz, double yy, double guessXX) const = 0;
        virtual double calcInverseY(double xx, double zz, double guessYY) const = 0;

    protected:
        std::vector<int> gridSizeX;
        std::vector<int> gridSizeY;
        std::vector<double> boundaryX;
        std::vector<double> boundaryY;
        std::vector<double> transformDx;
        std::vector<double> transformDy;
        Transformation transformationX;
        Transformation transformationY;

        std::vector<double> createInterpolationAxis(std::vector<int> gridSize, std::vector<double> boundary, Transformation transformation) const;

        double transform(double x, Transformation transformation) const;
        double backTransform(double x, Transformation transformation) const;

        
};

#endif // INTERPOLATION_HPP

std::vector<double> Interpolation::createInterpolationAxis(std::vector<int> gridSize, std::vector<double> boundary, Transformation transformation) const
{
    for (int i = 0; i < boundary.size(); i++)
    {
        boundary[i] = transform(boundary[i], transformation);
    }

    if(boundary[0] > boundary[1])
    {
        std::reverse(boundary.begin(), boundary.end());
    }

    double size = 0.0;
    for (int i = 0; i < gridSize.size(); i++)
    {
        size += gridSize[i];
    }

    std::vector<double> out(size+1);

    std::vector<double> dx = std::vector<double>(gridSize.size());

    for (int i = 0; i < dx.size(); i++)
    {
        dx[i] = (boundary[i+1] - boundary[i])/gridSize[i];
    }

    int idx = 0;
    out[idx] = backTransform(boundary[0], transformation);
    idx++;

    for (int ii = 0; ii < gridSize.size(); ii++)
    {
        for (int i = 0; i < gridSize[ii]; i++)
        {
            out[idx] = backTransform(boundary[ii] + dx[ii]*(i+1), transformation);
            idx++;
        }        
    }

    if(out[0] > out[1])
    {
        std::reverse(out.begin(), out.end());
    }

    return out;
}

double Interpolation::transform(double x, Transformation transformation) const
{
    switch(transformation) 
    {
        case LOG:
            return log(x);
        case LOG10:
            return log10(x);
        case LOGINV:
            return log(1/x);
    }

    return x;
}

double Interpolation::backTransform(double x, Transformation transformation) const
{
    switch(transformation) 
    {
        case LOG:
            return exp(x);
        case LOG10:
            return pow(x, 10.0);
        case LOGINV:
            return 1.0/exp(x);
    }

    return x;
}