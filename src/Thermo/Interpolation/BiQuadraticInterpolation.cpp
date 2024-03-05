#include <iostream>
#include "BiQuadraticInterpolation.hpp"

BiQuadraticInterpolation::BiQuadraticInterpolation(std::vector<double> xx,
                                                   std::vector<double> yy,
                                                   std::function<double(double, double)> f) : n(xx.size()), m(yy.size()), x(xx), y(yy), coeffs(n*m), Interpolation()
{
    dx = std::vector<double>(n+1);
    dy = std::vector<double>(m+1);
    z = std::vector<double>(n+1);
    t = std::vector<double>(m+1);

    for (int i = 1; i < dx.size()-1; i++)
    {
        dx[i] = x[i] - x[i-1];
    }
    dx[0] = dx[1];
    dx[n] = dx[n-1];

    for (int j = 1; j < dy.size()-1; j++)
    {
        dy[j] = y[j] - y[j-1];
    }
    dy[0] = dy[1];
    dy[m] = dy[m-1];

    z[0] = x[0] - 0.5*dx[0];
    for (int i = 0; i < z.size()-1; i++)
    {
        z[i+1] = x[i] + 0.5*dx[i+1];
    }

    t[0] = y[0] - 0.5*dy[0];
    for (int j = 0; j < t.size()-1; j++)
    {
        t[j+1] = y[j] + 0.5*dy[j+1];
    }

    calcCoeffs(f);
}



BiQuadraticInterpolation::BiQuadraticInterpolation(std::vector<double> xx,
                                                   std::vector<double> yy,
                                                   std::function<double(double, double)> f,
                                                   std::function<double(double, double)> fx,
                                                   std::function<double(double, double)> fy,
                                                   std::function<double(double, double)> fxy) : n(xx.size()), m(yy.size()), x(xx), y(yy), coeffs(n*m), Interpolation()
{
    dx = std::vector<double>(n+1);
    dy = std::vector<double>(m+1);
    z = std::vector<double>(n+1);
    t = std::vector<double>(m+1);

    for (int i = 1; i < dx.size()-1; i++)
    {
        dx[i] = x[i] - x[i-1];
    }
    dx[0] = dx[1];
    dx[n] = dx[n-1];

    for (int j = 1; j < dy.size()-1; j++)
    {
        dy[j] = y[j] - y[j-1];
    }
    dy[0] = dy[1];
    dy[m] = dy[m-1];

    z[0] = x[0] - 0.5*dx[0];
    for (int i = 0; i < z.size()-1; i++)
    {
        z[i+1] = x[i] + 0.5*dx[i+1];
    }

    t[0] = y[0] - 0.5*dy[0];
    for (int j = 0; j < t.size()-1; j++)
    {
        t[j+1] = y[j] + 0.5*dy[j+1];
    }

    calcCoeffs(f, fx, fy, fxy);
}


void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f,
                                          std::function<double(double, double)> fx,
                                          std::function<double(double, double)> fy,
                                          std::function<double(double, double)> fxy)
{
    Matrix u(n, m);
    Matrix p(n+1, m);
    Matrix q(n, m+1);
    Matrix r(n+1, m+1);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            u[i][j] = f(x[i], y[j]);
        }        
    }

    for (int j = 0; j < m; j++)
    {
        p[0][j] = fx(z[0], y[j]);
        p[n][j] = fx(z[n], y[j]);
    }

    for (int i = 0; i < n; i++)
    {
        q[i][0] = fy(x[i], t[0]);
        q[i][m] = fy(x[i], t[m]);
    }

    r[0][0] = fxy(z[0], t[0]);
    r[n][0] = fxy(z[n], t[0]);
    r[0][m] = fxy(z[0], t[m]);
    r[n][m] = fxy(z[n], t[m]);
    

    /////// solution for p
    for (int j = 0; j < m; j++)
    {
        std::vector<double> U(n-2);
        std::vector<double> D(n-1);
        std::vector<double> L(n-2);
        std::vector<double> b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*p[0][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*p[n][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        std::vector<double> pj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < pj.size()+1; i++)
        {
            p[i][j] = pj[i-1];
        }        
    }

    ///////// solution for q
    for (int i = 0; i < n; i++)
    {
        std::vector<double> U(m-2);
        std::vector<double> D(m-1);
        std::vector<double> L(m-2);
        std::vector<double> b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*q[i][0] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*q[i][m] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        std::vector<double> qi = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < qi.size()+1; j++)
        {
            q[i][j] = qi[j-1];
        }        
    }

    
    ///////// solution for ri0 and rim
    for (int j = 0; j < m+1; j += m)
    {
        std::vector<double> U(n-2);
        std::vector<double> D(n-1);
        std::vector<double> L(n-2);
        std::vector<double> b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*r[i][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*r[i+2][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        std::vector<double> rj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < rj.size()+1; i++)
        {
            r[i][j] = rj[i-1];
        }
    }


    //////solution for r
    for (int i = 0; i < n+1; i++)
    {
        std::vector<double> U(m-2);
        std::vector<double> D(m-1);
        std::vector<double> L(m-2);
        std::vector<double> b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*r[i][0] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*r[i][m] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        std::vector<double> ri = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < ri.size()+1; j++)
        {
            r[i][j] = ri[j-1];
        }        
    }

    int idx = 0;
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            Mat3x3 VxInv({0 , 1.0, 0,
                          dx[i+1]/(dx[i] + dx[i+1]), 0, dx[i]/(dx[i] + dx[i+1]),
                          -1.0/(dx[i] + dx[i+1]), 0, 1.0/(dx[i] + dx[i+1])});

            Mat3x3 VyInv({0 , 1.0, 0,
                          dy[j+1]/(dy[j] + dy[j+1]), 0, dy[j]/(dy[j] + dy[j+1]),
                          -1.0/(dy[j] + dy[j+1]), 0, 1.0/(dy[j] + dy[j+1])});

            Mat3x3 C({r[i][j], p[i][j], r[i][j+1],
                      q[i][j], u[i][j], q[i][j+1],
                      r[i+1][j], p[i+1][j], r[i+1][j+1]});

            coeffs[idx] = VyInv*C*(VyInv.transpose());
            idx++;
        }
    }
}


////////////////////

void BiQuadraticInterpolation::calcCoeffs(std::function<double(double, double)> f)
{
    Matrix u(n, m);
    Matrix p(n+1, m);
    Matrix q(n, m+1);
    Matrix r(n+1, m+1);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            u[i][j] = f(x[i], y[j]);
        }        
    }

    const double divider = 100.0;

    for (int j = 0; j < m; j++)
    {
        p[0][j] = calcDerivativeX(f, z[0], y[j], dx[0]/divider);
        p[n][j] = calcDerivativeX(f, z[n], y[j], dx[n]/divider);
    }

    for (int i = 0; i < n; i++)
    {
        q[i][0] = calcDerivativeY(f, x[i], t[0], dy[0]/divider);
        q[i][m] = calcDerivativeY(f, x[i], t[m], dy[m]/divider);
    }

    r[0][0] = calcDerivativeXY(f, z[0], t[0], dx[0]/divider, dy[0]/divider);
    r[n][0] = calcDerivativeXY(f, z[n], t[0], dx[n]/divider, dy[0]/divider);
    r[0][m] = calcDerivativeXY(f, z[0], t[m], dx[0]/divider, dy[m]/divider);
    r[n][m] = calcDerivativeXY(f, z[n], t[m], dx[n]/divider, dy[m]/divider);
    

    /////// solution for p
    for (int j = 0; j < m; j++)
    {
        std::vector<double> U(n-2);
        std::vector<double> D(n-1);
        std::vector<double> L(n-2);
        std::vector<double> b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*p[0][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*p[n][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        std::vector<double> pj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < pj.size()+1; i++)
        {
            p[i][j] = pj[i-1];
        }        
    }

    ///////// solution for q
    for (int i = 0; i < n; i++)
    {
        std::vector<double> U(m-2);
        std::vector<double> D(m-1);
        std::vector<double> L(m-2);
        std::vector<double> b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*q[i][0] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*q[i][m] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        std::vector<double> qi = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < qi.size()+1; j++)
        {
            q[i][j] = qi[j-1];
        }        
    }

    
    ///////// solution for ri0 and rim
    for (int j = 0; j < m+1; j += m)
    {
        std::vector<double> U(n-2);
        std::vector<double> D(n-1);
        std::vector<double> L(n-2);
        std::vector<double> b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*r[i][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*r[i+2][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        std::vector<double> rj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < rj.size()+1; i++)
        {
            r[i][j] = rj[i-1];
        }
    }


    //////solution for r
    for (int i = 0; i < n+1; i++)
    {
        std::vector<double> U(m-2);
        std::vector<double> D(m-1);
        std::vector<double> L(m-2);
        std::vector<double> b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*r[i][0] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*r[i][m] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        std::vector<double> ri = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < ri.size()+1; j++)
        {
            r[i][j] = ri[j-1];
        }        
    }

    int idx = 0;
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            Mat3x3 VxInv({0 , 1.0, 0,
                          dx[i+1]/(dx[i] + dx[i+1]), 0, dx[i]/(dx[i] + dx[i+1]),
                          -1.0/(dx[i] + dx[i+1]), 0, 1.0/(dx[i] + dx[i+1])});

            Mat3x3 VyInv({0 , 1.0, 0,
                          dy[j+1]/(dy[j] + dy[j+1]), 0, dy[j]/(dy[j] + dy[j+1]),
                          -1.0/(dy[j] + dy[j+1]), 0, 1.0/(dy[j] + dy[j+1])});

            Mat3x3 C({r[i][j], p[i][j], r[i][j+1],
                      q[i][j], u[i][j], q[i][j+1],
                      r[i+1][j], p[i+1][j], r[i+1][j+1]});

            coeffs[idx] = VyInv*C*(VyInv.transpose());
            idx++;
        }
    }
}

//////////////////

double BiQuadraticInterpolation::calc(double xx, double yy) const
{
    int i = 0;
    int found = 2;
    while(i < z.size()-1)
    {
        if(xx >= z[i] && xx <= z[i+1])
        {
            found--;
            break;
        }
        i = i + 1;
    }

    int j = 0;
    while(t.size()-1)
    {
        if(yy >= t[j] && yy <= t[j+1])
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

        Mat3x3 coeff = coeffs[j*n + i];

        return coeff[0][0] + w*(coeff[0][1] + w*coeff[0][2]) + v*(coeff[1][0] + w*(coeff[1][1] + w*coeff[1][2]) + v*(coeff[2][0] + w*(coeff[2][1] + w*coeff[2][2])));
    }

    std::cout << "ERROR: Interpolation out of range" << std::endl;
    return 0.0;
}


std::vector<double> BiQuadraticInterpolation::solveTridiagonal(const std::vector<double>& L,
                                                               const std::vector<double>& D,
                                                               const std::vector<double>& U,
                                                               const std::vector<double>& b) const
{
    int n = b.size();
    std::vector<double> out(n);
                                                                                                                                                                       
    std::vector<double> UStar(n-1, 0.0);
    std::vector<double> bStar(n, 0.0);
                                                                                                                                                    
    UStar[0] = U[0] / D[0];
    bStar[0] = b[0] / D[0];
                                                                                                                                            
    for (int i = 1; i < n-1; i++)
    {
        double m = 1.0/(D[i] - L[i-1]*UStar[i-1]);
        UStar[i] = U[i] * m;
        bStar[i] = (b[i] - L[i-1]*bStar[i-1])*m;
    }
    bStar[n-1] = (b[n-1] - L[n-2]*bStar[n-2])*(1.0/(D[n-1] - L[n-2] * UStar[n-2]));

    out[n-1] = bStar[n-1];
                                                                                                                                            
    for (int i = n-2; i >= 0; i--)
    {
        out[i] = bStar[i] - UStar[i]*out[i+1];
    }

    return out;
}

double BiQuadraticInterpolation::calcDerivativeX(std::function<double(double, double)> f, double x, double y, double step) const
{
    return (f(x+step, y) - f(x-step, y))/(2.0*step);
}

double BiQuadraticInterpolation::calcDerivativeY(std::function<double(double, double)> f, double x, double y, double step) const
{
    return (f(x, y+step) - f(x, y-step))/(2.0*step);
}

double BiQuadraticInterpolation::calcDerivativeXY(std::function<double(double, double)> f, double x, double y, double stepX, double stepY) const
{
    return (f(x+stepX, y+stepY) - f(x+stepX, y-stepY) - f(x-stepX, y+stepY) + f(x-stepX, y-stepY))/(4.0*stepX*stepY);
}