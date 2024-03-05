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


double BiLinearInterpolation::calc(double xx, double yy) const
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
        double v = xx - x[i];
        double w = yy - y[j];

        std::array<double, 4> nodeCoeffs = coeffs[j*n + i];

        return nodeCoeffs[0] + nodeCoeffs[1]*w + nodeCoeffs[2]*v + nodeCoeffs[3]*v*w;
    }

    std::cout << "ERROR: Interpolation out of range; " << xx << ", " << yy << std::endl;
    return 0.0;

}