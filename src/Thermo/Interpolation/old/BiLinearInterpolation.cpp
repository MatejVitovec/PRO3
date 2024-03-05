#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "BiLinearInterpolation.hpp"


void BiLinearInterpolation::calcCoeffs(std::function<double(double, double)> f)
{
    std::vector<double> knotValues = std::vector<double>((sizeX+1)*(sizeY+1));

    int cnt = 0;

    int jStart = 0;
    int jEnd = ySizes[0] + 1;
    int jIndex = 0;
    for (int jj = 0; jj < ySizes.size(); jj++)
    {
        for (int j = jStart; j < jEnd; j++)
        {
            int iStart = 0;
            int iEnd = xSizes[0] + 1;
            int iIndex = 0;
            for (int ii = 0; ii < xSizes.size(); ii++)
            {
                for (int i = iStart; i < iEnd; i++)
                {
                    //std::pair<double, double> xyBackTransformed = backTransform(xBoundary[ii] + dx[ii]*iIndex, yBoundary[jj] + dy[jj]*jIndex);

                    knotValues[j*(sizeY+1) + i] = f(backTransformX(xBoundary[ii] + dx[ii]*iIndex), backTransformY(yBoundary[jj] + dy[jj]*jIndex));
                    cnt++;
                    iIndex++;
                }
                iStart = iEnd;
                iEnd = iStart + xSizes[ii+1];
                iIndex = 1;
            }
            jIndex++;
        }
        jStart = jEnd;
        jEnd = jStart + ySizes[jj+1];
        jIndex = 1;

        std::cout << cnt << std::endl;
    }


    int firstCoeffIndex = 0;
    for (int j = 0; j < sizeY; j++)
    {
        for (int i = 0; i < sizeX; i++)
        {
            coeffs[firstCoeffIndex]   = knotValues[j*(sizeY+1) + i];
            coeffs[firstCoeffIndex+1] = knotValues[j*(sizeY+1) + (i+1)] - knotValues[j*(sizeY+1) + i];
            coeffs[firstCoeffIndex+2] = knotValues[(j+1)*(sizeY+1) + i] - knotValues[j*(sizeY+1) + i];
            coeffs[firstCoeffIndex+3] = knotValues[(j+1)*(sizeY+1) + (i+1)] - knotValues[(j+1)*(sizeY+1) + i] - knotValues[j*(sizeY+1) + (i+1)] + knotValues[j*(sizeY+1) + i];
            firstCoeffIndex += 4;
        }
    }
}

std::array<double, 4> BiLinearInterpolation::getNodeCoeffs(int i, int j) const
{
    std::array<double, 4> out;
    std::copy_n(coeffs.begin() + j*(sizeX*4) + i*4, 4, out.begin());
    return out;    
}

std::tuple<double, double, double, double> BiLinearInterpolation::getNodeXYdXdY(int i, int j) const
{
    int ii;
    int jj;

    int nXAux = 0;
    for (ii = 0; ii < xSizes.size(); ii++)
    {
        if (xSizes[ii] + nXAux > i) break;
        nXAux += xSizes[ii];
    }

    int nYAux = 0;
    for (jj = 0; jj < ySizes.size(); jj++)
    {
        if (ySizes[jj] + nYAux > j) break;
        nYAux += ySizes[jj];
    }
    
    return {xBoundary[ii] + dx[ii]*(i - nXAux), yBoundary[jj] + dy[jj]*(j - nYAux), dx[ii], dy[jj]};
}


double BiLinearInterpolation::calc(double x, double y) const
{
    x = transformX(x);
    y = transformY(y);
    
    if(!checkInterval(x, y))
    {
        std::cout << "Error: Interpolation interval exceeded" << std::endl;
        return 0.0;
    }

    std::pair<int, int> position = findPosition(x, y);
    int i = position.first;
    int j = position.second;

    std::tuple<double, double, double, double> xydxdy = getNodeXYdXdY(i, j);
    double xi = std::get<0>(xydxdy);
    double yj = std::get<1>(xydxdy);
    double dxi = std::get<2>(xydxdy);
    double dyj = std::get<3>(xydxdy);

    std::array<double, 4>  nodeCoeffs = getNodeCoeffs(i, j);

    return nodeCoeffs[0] + nodeCoeffs[1]*((x-xi)/dxi) + nodeCoeffs[2]*((y-yj)/dyj) + nodeCoeffs[3]*(((x-xi)*(y-yj))/(dxi*dyj));
}

double BiLinearInterpolation::calcInverseX(double y, double z, double xGuess) const
{
    xGuess = transformX(xGuess);
    y = transformY(y);

    std::pair<int, int> position = findPosition(xGuess, y);
    int i = position.first;
    int j = position.second;
    
    double xTransformed;

    while (true)
    {
        std::tuple<double, double, double, double> xydxdy = getNodeXYdXdY(i, j);
        double xi = std::get<0>(xydxdy);
        double yj = std::get<1>(xydxdy);
        double dxi = std::get<2>(xydxdy);
        double dyj = std::get<3>(xydxdy);

        std::array<double, 4>  nodeCoeffs = getNodeCoeffs(i, j);

        xTransformed = (z - nodeCoeffs[0] + nodeCoeffs[1]*(xi/dxi) - nodeCoeffs[2]*((y - yj)/dyj) + nodeCoeffs[3]*xi*((y - yj)/(dxi*dyj)))/
                       ((nodeCoeffs[1] + nodeCoeffs[3]*((y-yj)/dyj))/dxi);

        if (xi < xi+dxi)
        {
            if (xTransformed < xi) { i--; }
            else if(xTransformed > xi+dxi) { i++; }
            else { break; }
        }
        else
        {
            if (xTransformed > xi) { i--; }
            else if(xTransformed < xi+dxi) { i++; }
            else { break; }
        }
    }

    return backTransformX(xTransformed);
}

double BiLinearInterpolation::calcInverseY(double x, double z, double yGuess) const
{
    x = transformX(x);
    yGuess = transformY(yGuess);

    std::pair<double, double> position = findPosition(x, yGuess);
    int i = position.first;
    int j = position.second;
    
    double yTransformed;

    while (true)
    {
        std::tuple<double, double, double, double> xydxdy = getNodeXYdXdY(i, j);
        double xi = std::get<0>(xydxdy);
        double yj = std::get<1>(xydxdy);
        double dxi = std::get<2>(xydxdy);
        double dyj = std::get<3>(xydxdy);

        std::array<double, 4>  nodeCoeffs = getNodeCoeffs(i, j);

        yTransformed = (z - nodeCoeffs[0] + nodeCoeffs[1]*((x-xi)/dxi) - nodeCoeffs[2]*((yj)/dyj) + nodeCoeffs[3]*xi*((x - xi)/(dxi*dyj)))/
                       ((nodeCoeffs[2] + nodeCoeffs[3]*((x-xi)/dxi))/dyj);

        if (yj < yj+dyj)
        {
            if (yTransformed < yj) { j--; }
            else if(yTransformed > yj+dyj) { j++; }
            else { break; }
        }
        else
        {
            if (yTransformed > yj) { j--; }
            else if(yTransformed < yj+dyj) { j++; }
            else { break; }
        }
    }

    return backTransformY(yTransformed);
}