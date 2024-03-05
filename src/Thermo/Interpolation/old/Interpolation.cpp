#include <cmath>
#include <iostream>
#include <algorithm>

#include "Interpolation.hpp"

Interpolation::Interpolation(std::vector<int> xSizes_,
                             std::vector<int> ySizes_,
                             std::vector<double> xBoundary_,
                             std::vector<double> yBoundary_,
                             TransformationType xTransformation_,
                             TransformationType yTransformation_) : xSizes(xSizes_),
                                                                    ySizes(ySizes_),
                                                                    xTransformation(xTransformation_), 
                                                                    yTransformation(yTransformation_)
{
    xBoundary = xBoundary_;
    yBoundary = yBoundary_;

    for (int i = 0; i < xBoundary_.size(); i++)
    {
        switch(xTransformation) 
        {
            case LOG:
                xBoundary[i] = log(xBoundary_[i]);
                break;
            case LOG10:
                xBoundary[i] = log10(xBoundary_[i]);
                break;
            case LOGINV:
                xBoundary[i] = log(1/xBoundary_[i]);
                break;
        }
    }

    for (int i = 0; i < yBoundary_.size(); i++)
    {
        switch(yTransformation) 
        {
            case LOG:
                yBoundary[i] = log(yBoundary_[i]);
                break;
            case LOG10:
                yBoundary[i] = log10(yBoundary_[i]);
                break;
            case LOGINV:
                yBoundary[i] = log(1/yBoundary_[i]);
                break;
        }
    }

    if(xBoundary[0] > xBoundary[1])
    {
        std::reverse(xBoundary.begin(), xBoundary.end());
    }
    if(yBoundary[0] > yBoundary[1])
    {
        std::reverse(yBoundary.begin(), yBoundary.end());
    }


    sizeX = 0.0;
    for (int i = 0; i < xSizes.size(); i++)
    {
        sizeX += xSizes[i];
    }

    sizeY = 0.0;
    for (int j = 0; j < ySizes.size(); j++)
    {
        sizeY += ySizes[j];
    }

    dx = std::vector<double>(xSizes.size());
    dy = std::vector<double>(ySizes.size());

    for (int i = 0; i < dx.size(); i++)
    {
        dx[i] = (xBoundary[i+1] - xBoundary[i])/xSizes[i];
    }

    for (int i = 0; i < dy.size(); i++)
    {
        dy[i] = (yBoundary[i+1] - yBoundary[i])/ySizes[i];
    }
}

double Interpolation::xMax() const
{
    return xBoundary.back();
}

double Interpolation::yMax() const
{
    return yBoundary.back();
}

double Interpolation::transformX(double x) const
{
    switch(xTransformation) 
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

double Interpolation::transformY(double y) const
{
    switch(yTransformation) 
    {
        case LOG:
            return log(y);
        case LOG10:
            return log10(y);
        case LOGINV:
            return log(1/y);
    }

    return y;
}

double Interpolation::backTransformX(double x) const
{
    switch(xTransformation) 
    {
        case LOG:
            return exp(x);
        case LOG10:
            return pow(x, 10.0);
        case LOGINV:
            return 1/exp(x);
    }

    return x;
}

double Interpolation::backTransformY(double y) const
{
    switch(yTransformation) 
    {
        case LOG:
            return exp(y);
        case LOG10:
            return pow(y, 10.0);
        case LOGINV:
            return 1/exp(y);
    }

    return y;
}

bool Interpolation::checkInterval(double x, double y) const
{
    if (x < xBoundary[0] || x > xMax() || y < yBoundary[0] || y > yMax()) //TODO obousmerny interval - transformace muze otocit znamenka
    {
        return false;
    }
    
    return true;
}

std::pair<int, int> Interpolation::findPosition(double x, double y) const
{
    int ii;
    int jj;
    for (ii = 0; ii < xSizes.size(); ii++) { if (x < xBoundary[ii+1]) break; }
    for (jj = 0; jj < ySizes.size(); jj++) { if (y < yBoundary[jj+1]) break; }

    int xIdx = std::floor((x - xBoundary[ii])/dx[ii]);
    int yIdx = std::floor((y - yBoundary[jj])/dy[jj]);

    return std::make_pair(xIdx + ((ii == 0) ? 0 : xSizes[ii-1]),
                          yIdx + ((jj == 0) ? 0 : ySizes[jj-1]));
}


