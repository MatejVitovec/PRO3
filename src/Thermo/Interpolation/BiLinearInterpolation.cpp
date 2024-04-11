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
    std::vector<double> boundaryXAux = boundaryX;
    std::vector<double> boundaryYAux = boundaryY;
    for (int i = 0; i < boundaryXAux.size(); i++)
    {
        boundaryXAux[i] = transform(boundaryXAux[i], transformationX);
    }

    for (int j = 0; j < boundaryYAux.size(); j++)
    {
        boundaryYAux[j] = transform(boundaryYAux[j], transformationY);
    }

    if(boundaryXAux[0] > boundaryXAux[1])
    {
        std::reverse(boundaryXAux.begin(), boundaryXAux.end());
    }

    if(boundaryYAux[0] > boundaryYAux[1])
    {
        std::reverse(boundaryYAux.begin(), boundaryYAux.end());
    }

    double sizeX = 0.0;
    for (int i = 0; i < gridSizeX.size(); i++)
    {
        sizeX += gridSizeX[i];
    }

    double sizeY = 0.0;
    for (int j = 0; j < gridSizeY.size(); j++)
    {
        sizeY += gridSizeY[j];
    }

    x = std::vector<double>(sizeX+1);
    y = std::vector<double>(sizeY+1);

    dxTransf = std::vector<double>(gridSizeX.size());
    dyTransf = std::vector<double>(gridSizeY.size());

    for (int i = 0; i < dxTransf.size(); i++)
    {
        dxTransf[i] = (boundaryXAux[i+1] - boundaryXAux[i])/gridSizeX[i];
    }

    for (int j = 0; j < dyTransf.size(); j++)
    {
        dyTransf[j] = (boundaryYAux[j+1] - boundaryYAux[j])/gridSizeY[j];
    }

    int idx = 0;
    x[idx] = backTransform(boundaryXAux[0], transformationX);
    idx++;

    for (int ii = 0; ii < gridSizeX.size(); ii++)
    {
        for (int i = 0; i < gridSizeX[ii]; i++)
        {
            x[idx] = backTransform(boundaryXAux[ii] + dxTransf[ii]*(i+1), transformationX);
            idx++;
        }        
    }

    idx = 0;
    y[idx] = backTransform(boundaryYAux[0], transformationY);
    idx++;

    for (int jj = 0; jj < gridSizeY.size(); jj++)
    {
        for (int j = 0; j < gridSizeY[jj]; j++)
        {
            y[idx] = backTransform(boundaryYAux[jj] + dyTransf[jj]*(j+1), transformationY);
            idx++;
        }        
    }

    if(x[0] > x[1]) { std::reverse(x.begin(), x.end()); }
    if(y[0] > y[1]) { std::reverse(y.begin(), y.end()); }


    n = x.size()-1;
    m = y.size()-1;

    dx = std::vector<double>(n);
    dy = std::vector<double>(m);

    coeffs = std::vector<std::array<double, 4>>(m*n);

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

double BiLinearInterpolation::calcFastFind(double xx, double yy) const
{
    std::pair<int, int> position = fastFindPosition(xx, yy);

    double v = xx - x[position.first];
    double w = yy - y[position.second];

    std::array<double, 4> nodeCoeffs = coeffs[position.second*n + position.first];

    return nodeCoeffs[0] + nodeCoeffs[1]*w + nodeCoeffs[2]*v + nodeCoeffs[3]*v*w;
}

double BiLinearInterpolation::calcFastFindDiffX(double xx, double yy) const
{
    std::pair<int, int> position = findPosition(xx, yy);

    double v = xx - x[position.first];
    double w = yy - y[position.second];

    std::array<double, 4> nodeCoeffs = coeffs[position.second*n + position.first];

    return nodeCoeffs[2] + nodeCoeffs[3]*w;
}

double BiLinearInterpolation::calcFastFindDiffY(double xx, double yy) const
{
    std::pair<int, int> position = findPosition(xx, yy);

    double v = xx - x[position.first];
    double w = yy - y[position.second];

    std::array<double, 4> nodeCoeffs = coeffs[position.second*n + position.first];

    return nodeCoeffs[1] + nodeCoeffs[3]*v;
}

double BiLinearInterpolation::calcInverseX(double zz, double yy, double guessXX) const
{
    std::pair<int, int> position = findPosition(guessXX, yy);

    std::array<double, 4> coeff = coeffs[position.second*n + position.first];

    double w = yy - y[position.first];

    return (zz - coeff[0] - coeff[1]*w)/(coeff[2] + coeff[3]*w) + x[position.second];
}

double BiLinearInterpolation::calcInverseY(double xx, double zz, double guessYY) const
{
    std::pair<int, int> position = findPosition(xx, guessYY);

    std::array<double, 4> coeff = coeffs[position.second*n + position.first];

    double v = xx - x[position.first];

    return (zz - coeff[0] - coeff[2]*v)/(coeff[1] + coeff[3]*v) + y[position.second];
}