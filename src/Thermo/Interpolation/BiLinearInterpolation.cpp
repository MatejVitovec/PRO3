#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "BiLinearInterpolation.hpp"

BiLinearInterpolation::BiLinearInterpolation(std::vector<double> xx,
                                             std::vector<double> yy,
                                             std::function<double(double, double)> f) : n(xx.size()-1), m(yy.size()-1), x(xx), y(yy), coeffs(n*m), Interpolation()
{
    dx = std::vector<double>(n);
    dy = std::vector<double>(m);

    for (int i = 0; i < dx.size(); i++)
    {
        dx[i] = x[i+1] - x[i];
    }


    for (int j = 0; j < dy.size(); j++)
    {
        dy[j] = y[j+1] - y[j];
    }

    calcCoeffs(f);
}

BiLinearInterpolation::BiLinearInterpolation(std::vector<int> gridSizeX_, std::vector<int> gridSizeY_,
                                             std::vector<double> boundaryX_, std::vector<double> boundaryY_,
                                             Transformation transformationX_, Transformation transformationY_,
                                             std::function<double(double, double)> f) : n(0), m(0), x(0), y(0), coeffs(0),
                                                Interpolation(gridSizeX_, gridSizeY_, boundaryX_, boundaryY_, transformationX_, transformationY_)
{
    //TODO
    dx = std::vector<double>(n);
    dy = std::vector<double>(m);

    for (int i = 0; i < dx.size(); i++)
    {
        dx[i] = x[i+1] - x[i];
    }


    for (int j = 0; j < dy.size(); j++)
    {
        dy[j] = y[j+1] - y[j];
    }

    calcCoeffs(f);
}


void BiLinearInterpolation::calcCoeffs(std::function<double(double, double)> f)
{
    int idx = 0;
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            coeffs[idx] = std::array<double, 4>({f(x[i], y[j]),
                                                 (f(x[i], y[j+1]) - f(x[i], y[j]))/dy[j],
                                                 (f(x[i+1], y[j]) - f(x[i], y[j]))/dx[i],
                                                 (f(x[i+1], y[j+1]) - f(x[i], y[j+1]) - f(x[i+1], y[j]) + f(x[i], y[j]))/(dx[i]*dy[j])});
            idx++;
        }        
    }
}

std::pair<int, int> BiLinearInterpolation::findPosition(double xx, double yy) const
{
    int i = 0;
    int found = 2;
    while(i < x.size()-1)
    {
        if(xx >= x[i] && xx <= x[i+1])
        {
            found--;
            break;
        }
        i = i + 1;
    }

    int j = 0;
    while(y.size()-1)
    {
        if(yy >= y[j] && yy <= y[j+1])
        {
            found--;
            break;
        }
        j = j + 1;
    }

    if(found == 0)
    {
        return std::pair<int, int>(i, j);
    }

    std::cout << "ERROR: Interpolation out of range" << std::endl;
    return std::pair<int, int>(0, 0);
}

std::pair<int, int> BiLinearInterpolation::fastFindPosition(double xx, double yy) const
{
    int shiftIdxX = 0;
    int shiftIdxY = 0;

    int ii;
    for (ii = 0; ii < gridSizeX.size(); ii++)
    {
        if (xx < boundaryX[ii+1]) { break; }
        shiftIdxX += gridSizeX[ii];
    }
    int jj;
    for (jj = 0; jj < gridSizeY.size(); jj++)
    {
        if (yy < boundaryY[jj+1]) { break; }
        shiftIdxY += gridSizeY[jj];
    }
    
    return std::pair<int, int>(std::floor((abs(transform(xx, transformationX) - transform(boundaryX[ii], transformationX)))/dxTransf[ii]) + shiftIdxX,
                               std::floor((abs(transform(yy, transformationY) - transform(boundaryY[jj], transformationY)))/dyTransf[jj]) + shiftIdxY);
}

double BiLinearInterpolation::calc(double xx, double yy) const
{
    std::pair<int, int> position = findPosition(xx, yy);

    double v = xx - x[position.first];
    double w = yy - y[position.second];

    std::array<double, 4> nodeCoeffs = coeffs[position.second*n + position.first];

    return nodeCoeffs[0] + nodeCoeffs[1]*w + nodeCoeffs[2]*v + nodeCoeffs[3]*v*w;
}

double BiLinearInterpolation::calcInverseX(double zz, double yy, double guessXX) const
{
    //TODO
    return 0.0;
}

double BiLinearInterpolation::calcInverseY(double xx, double zz, double guessYY) const
{
    //TODO
    return 0.0;
}