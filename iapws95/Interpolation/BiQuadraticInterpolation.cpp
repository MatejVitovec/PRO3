#include <cmath>
#include <iostream>
#include <fstream>

#include "BiQuadraticInterpolation.hpp"

void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f)
{



    //TODO
}

void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f,
                                          std::function<double(double, double)> fx,
                                          std::function<double(double, double)> fy,
                                          std::function<double(double, double)> fxy)
{
    std::vector<double> nodeValues = std::vector<double>((sizeX)*(sizeY));
    std::vector<double> xDerivative0 = std::vector<double>(sizeX+1);
    std::vector<double> xDerivativeN = std::vector<double>(sizeX+1);
    std::vector<double> yDerivative0 = std::vector<double>(sizeY+1);
    std::vector<double> yDerivativeM = std::vector<double>(sizeY+1);
    std::array<double, 4> xyCornerDerivative;

    std::vector<double> xKnotVec(sizeX+1);
    std::vector<double> dxVec(sizeX);
    std::vector<double> yKnotVec(sizeY+1);
    std::vector<double> dyVec(sizeY);

    //calculate dx dy between knots 
    {
        int iStart = 0;
        int iEnd = xSizes[0] + 1;
        int iIndex = 0;
        for (int ii = 0; ii < xSizes.size(); ii++)
        {
            for (int i = iStart; i < iEnd; i++)
            {
                xKnotVec[i] = xBoundary[ii] + dx[ii]*iIndex;
                
                if(i > 0)
                {
                    dxVec[i-1] = dx[ii];                    
                }
                iIndex++;
            }
            iStart = iEnd;
            iEnd = iStart + xSizes[ii];
            iIndex = 0;
        }        

        int jStart = 0;
        int jEnd = ySizes[0] + 1;
        int jIndex = 0;
        for (int jj = 0; jj < ySizes.size(); jj++)
        {
            for (int j = jStart; j < jEnd; j++)
            {
                yKnotVec[j] = yBoundary[jj] + dy[jj]*jIndex;
                
                if(j > 0)
                {
                    dyVec[j-1] = dy[jj];
                }
                jIndex++;
            }
            jStart = jEnd;
            jEnd = jStart + ySizes[jj];
            jIndex = 0;
        }
    }

    //calculate nodes values
    for (int j = 0; j < sizeY; j++)
    {
        for (int i = 0; i < sizeX; i++)
        {
            nodeValues[j*sizeY + i] = f(backTransformX(xKnotVec[i] + dxVec[i]/2.0), backTransformY(yKnotVec[i] + dyVec[i]/2.0));
        }
    }

    //transform dx dy between knots to dx dy between nodes
    for (int i = 0; i < dxVec.size()-1; i++)
    {
        dxVec[i] = (dxVec[i] + dxVec[i+1])/2.0;
    }
    dxVec.pop_back();

    for (int j = 0; j < dyVec.size()-1; j++)
    {
        dyVec[j] = (dyVec[j] + dyVec[j+1])/2.0;
    }
    dyVec.pop_back();

    std::vector<double> p((sizeX+1)*(sizeY+1));
    std::vector<double> q((sizeX+1)*(sizeY+1));
    std::vector<double> r((sizeX+1)*(sizeY+1));

    for (int j = 0; j < sizeY+1; j++)
    {
        p[j*(sizeX+1) + 0] = fx(backTransformX(xKnotVec[0]), backTransformY(yKnotVec[j]));
        p[j*(sizeX+1) + sizeX] = fx(backTransformX(xKnotVec[xKnotVec.size()-1]), backTransformY(yKnotVec[j]));

        //xDerivative0[j] = fx(backTransformX(xKnotVec[0]), backTransformY(yKnotVec[j]));
        //xDerivativeN[j] = fx(backTransformX(xKnotVec[xKnotVec.size()-1]), backTransformY(yKnotVec[j]));
    }

    for (int i = 0; i < sizeX+1; i++)
    {
        q[0*(sizeX+1) + i] = fy(backTransformX(xKnotVec[i]), backTransformY(yKnotVec[0]));
        q[sizeY*(sizeX+1) + i] = fy(backTransformX(xKnotVec[i]), backTransformY(yKnotVec[yKnotVec.size()-1]));

        //yDerivative0[i] = fy(backTransformX(xKnotVec[i]), backTransformY(yKnotVec[0]));
        //yDerivativeM[i] = fy(backTransformX(xKnotVec[i]), backTransformY(yKnotVec[yKnotVec.size()-1]));
    }
    
    r[0*(sizeX+1) + 0] = fxy(backTransformX(xBoundary[0]), backTransformY(yBoundary[0]));
    r[0*(sizeX+1) + sizeX] = fxy(backTransformX(xBoundary[xSizes.size()]), backTransformY(yBoundary[0]));
    r[sizeY*(sizeX+1) + 0] = fxy(backTransformX(xBoundary[0]), backTransformY(yBoundary[ySizes.size()]));
    r[sizeY*(sizeX+1) + sizeX] = fxy(backTransformX(xBoundary[xSizes.size()]), backTransformY(yBoundary[ySizes.size()]));

    //xyCornerDerivative[0] = fxy(backTransformX(xBoundary[0]), backTransformY(yBoundary[0]));
    //xyCornerDerivative[1] = fxy(backTransformX(xBoundary[xSizes.size()]), backTransformY(yBoundary[0]));
    //xyCornerDerivative[2] = fxy(backTransformX(xBoundary[0]), backTransformY(yBoundary[ySizes.size()]));
    //xyCornerDerivative[3] = fxy(backTransformX(xBoundary[xSizes.size()]), backTransformY(yBoundary[ySizes.size()]));

    //create matrixes



    for (int j = 0; j < sizeY; j++)
    {
        std::vector<double> L(sizeX);
        std::vector<double> D(sizeX+1);
        std::vector<double> U(sizeX);
        std::vector<double> b(sizeX-1);

        for (int i = 1; i < sizeX; i++)
        {
            L[i-1] = 1/(dxVec[i] + dxVec[i+1]);
            D[i-1] = (2 + dxVec[i-1]/(dxVec[i-1] + dxVec[i]) + dxVec[i+1]/(dxVec[i] + dxVec[i+1]))/dxVec[i];
            U[i-1] = 1/(dxVec[i] + dxVec[i+1]);
            b[i-1] = 4*(nodeValues[j*sizeY + i+1] - nodeValues[j*sizeY + i])/(dxVec[i]*dxVec[i]);
        }
        D[sizeX-1] = (2 + dxVec[sizeX-1]/(dxVec[sizeX-1] + dxVec[sizeX]) + dxVec[sizeX+1]/(dxVec[sizeX] + dxVec[sizeX+1]))/dxVec[sizeX];
        //b[sizeX-2] = 4*(nodeValues[j*sizeY + i+1] - nodeValues[j*sizeY + i])/(dxVec[i]*dxVec[i]);
        

    }

    /*for (int i = 0; i < sizeX; i++)
    {
        int j = 0;
        rhsP[i] -= xDerivative0[j]/(dxVec[i-1] + dxVec[i]);
        tripletListP.push_back(T(i,   j, (2 + dxVec[i-1]/(dxVec[i-1] + dxVec[i]) + dxVec[i+1]/(dxVec[i] + dxVec[i+1]))));
        tripletListP.push_back(T(i+1, j, 1/(dxVec[i] + dxVec[i+1])));
        rhsP[i] += 4*(nodeValues[j*sizeY + i+1] - nodeValues[j*sizeY + i])/(dxVec[i]*dxVec[i]);
    }*/
    
    /////

    /*for (int j = 0; j < sizeY; j++)
    {
        for (int i = 1; i < sizeX-2; i++)
        {
            tripletListP.push_back(T(i-1, j, 1/(dxVec[i-1] + dxVec[i])));
            tripletListP.push_back(T(i,   j, (2 + dxVec[i-1]/(dxVec[i-1] + dxVec[i]) + dxVec[i+1]/(dxVec[i] + dxVec[i+1]))));
            tripletListP.push_back(T(i+1, j, 1/(dxVec[i] + dxVec[i+1])));
            rhsP[j] += 4*(nodeValues[j*sizeY + i+1] - nodeValues[j*sizeY + i])/(dxVec[i]*dxVec[i]);
        }
    }

    for (int j = 1; j < sizeY-2; j++)
    {
        for (int i = 0; i < sizeX; i++)
        {
            tripletListQ.push_back(T(i, j-1, 1/(dyVec[j-1] + dyVec[j])));
            tripletListQ.push_back(T(i, j,   (2 + dyVec[j-1]/(dyVec[j-1] + dyVec[j]) + dyVec[j+1]/(dyVec[j] + dyVec[j+1]))));
            tripletListQ.push_back(T(i, j+1, 1/(dyVec[j] + dyVec[j+1])));
            rhsQ[j] += 4*(nodeValues[(j+1)*sizeY + i] - nodeValues[j*sizeY + j])/(dyVec[j]*dyVec[j]);
        }
    }*/
    

    //solve matrixes


    //calculate coeffs

}

BiQuadraticInterpolation::BiQuadraticInterpolationCoeffs BiQuadraticInterpolation::getNodeCoeffs(int i, int j) const
{
    return coeffs[j*sizeX + i];
}


std::pair<double, double> BiQuadraticInterpolation::getNodeXY(int i, int j) const
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
    
    return {xBoundary[ii] + dx[ii]*(i - nXAux), yBoundary[jj] + dy[jj]*(j - nYAux)};
}


double BiQuadraticInterpolation::calc(double x, double y) const
{
    x = transformX(x);
    y = transformY(y);
    
    if(!checkInterval(x, y))
    {
        std::cout << "Error: Interpolation interval exceeded" << std::endl;
        return 0.0;
    }

    std::pair<double, double> position = findPosition(x, y);
    int i = position.first;
    int j = position.second;

    std::pair<double, double> xiyj = getNodeXY(i, j);
    double s = x - xiyj.first;
    double t = y - xiyj.second;

    BiQuadraticInterpolationCoeffs A = getNodeCoeffs(i, j);

    return A[0][0] + t*(A[0][1] + t*A[0][2]) + s*(A[1][0] + t*(A[1][1] + t*A[1][2]) + s*(A[2][0] + t*(A[2][1] + t*A[2][2])));
}

/*double BiQuadraticInterpolation::calcInverseX(double y, double z, double xGuess) const
{
    xGuess = transformX(xGuess);
    y = transformY(y);

    std::pair<double, double> position = findPosition(xGuess, y);
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

/*double BiQuadraticInterpolation::calcInverseY(double x, double z, double yGuess) const
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
}*/

std::vector<double> BiQuadraticInterpolation::solveTridiagonalSystem(const std::vector<double>& L, const std::vector<double>& D, const std::vector<double>& U, const std::vector<double>& b)
{
    int n = b.size();
    std::vector<double> out(n);
                                                                                                                                                                       
    std::vector<double> UStar(n, 0.0);
    std::vector<double> bStar(n, 0.0);
                                                                                                                                                    
    UStar[0] = U[0] / D[0];
    bStar[0] = b[0] / D[0];
                                                                                                                                            
    for (int i = 1; i < n; i++)
    {
        double m = 1.0 / (D[i] - L[i] * UStar[i-1]);
        UStar[i] = U[i] * m;
        bStar[i] = (b[i] - L[i] * bStar[i-1]) * m;
    }
                                                                                                                                            
    for (int i = n-1; i > 0; i--)
    {
        out[i] = bStar[i] - UStar[i] * b[i+1];
    }

    return out;
}