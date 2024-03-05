#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

#include <fenv.h>

#include "Iawps95.hpp"
#include "Interpolation/BiLinearInterpolation.hpp"


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    Iawps95 thermo = Iawps95();
    thermo.saturatedCurve();

    // [x, y] = [rho, e]

    std::vector<std::vector<double>> volumeT(200, std::vector<double>(200));
    std::vector<std::vector<double>> energy(200, std::vector<double>(200));

    std::vector<std::vector<double>> density(200, std::vector<double>(200));
    std::vector<std::vector<double>> pressure(200, std::vector<double>(200));
    std::vector<std::vector<double>> temperature(200, std::vector<double>(200));
    std::vector<std::vector<double>> soundSpeed(200, std::vector<double>(200));

    double dx = (log(0.008) - log(1188.87))/200;
    for (int i = 0; i < 200; i++)
    {
        double densityAux = log(1188.87) + i*dx;
        for (int j = 0; j < 200; j++)
        {
            volumeT[i][j] = densityAux;
            density[i][j] = 1/exp(volumeT[i][j]);
        }        
    }

    //dx = (log(0.008) - log(0.00102796))/200;
    //for (int i = 0; i < 200; i++)
    //{
    //    double densityAux = log(0.00102796) + i*dx;
    //    for (int j = 0; j < 200; j++)
    //    {
    //        volumeT[i+200][j] = densityAux;
    //        density[i+200][j] = 1/exp(volumeT[i+200][j]);
    //    }        
    //}

    double dy = (2650000.0 - 2005000.0)/100;
    for (int j = 0; j < 100; j++)
    {
        double energyAux = 2005000.0 + j*dy;
        for (int i = 0; i < 200; i++)
        {
            energy[i][j] = energyAux;
        }        
    }

    dy = (4085000.27 - 2650000.0)/100;
    for (int j = 0; j < 100; j++)
    {
        double energyAux = 2650000.0 + j*dy;
        for (int i = 0; i < 200; i++)
        {
            energy[i][j+100] = energyAux;
        }        
    }

    int cnt = 0;
    
    for (int i = 0; i < 200; i++)
    {
        for (int j = 0; j < 200; j++)
        {
            temperature[i][j] = thermo.implicitTFromRhoE(density[i][j], energy[i][j], (energy[i][j]*(1.32 - 1.0))/461.51805);
            pressure[i][j] = thermo.p(density[i][j], temperature[i][j]);
            soundSpeed[i][j] = sqrt(thermo.w2(density[i][j], temperature[i][j]));
            cnt++;
            std::cout << cnt << std::endl;
        }
    }

    std::ofstream f;
	f.open("interpData2.dat", std::ios::out);

    double maxTemp = 0;
    for (int j = 0; j < 200; j++)
    {
        for (int i = 0; i < 200; i++)
        {
            f << density[i][j] << " " << energy[i][j] << " " << temperature[i][j] << " " << pressure[i][j] << " " << soundSpeed[i][j] << "\n";
            maxTemp = std::max(temperature[i][j], maxTemp);
        }
    }

    f << std::endl;
    f.close();

    std::cout << maxTemp << std::endl;

    //Iawps95 thermo = Iawps95();
    
    //BiLinearInterpolation interp = BiLinearInterpolation();
    //BiLinearInterpolation interp = BiLinearInterpolation(std::vector<int>{200}, std::vector<int>{100, 100}, std::vector<double>{1/1188.87, 1/0.008}, std::vector<double>{2005000.0, 2650000.0, 4085000.27}, Interpolation::LOGINV, Interpolation::NONE);
    //std::unique_ptr<Interpolation> interp = std::make_unique<BiLinearInterpolation>(std::vector<int>{200}, std::vector<int>{100, 100}, std::vector<double>{1/1188.87, 0.008}, std::vector<double>{2005000.0, 2650000.0, 4085000.27}, Interpolation::LOGINV, Interpolation::NONE);
    //std::unique_ptr<Interpolation> interp = std::make_unique<BiLinearInterpolation>();

    /*interp.calcCoeffs([=](double density, double energy) { return thermo.implicitTFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

    double rho = 0.5;
    double energy = 3000000.0;

    double test = interp(rho, energy);
    double testIapws = thermo.implicitTFromRhoE(rho, energy, (energy*(1.32 - 1.0))/461.51805);

    double rhoInv = interp.calcInverseX(energy, test, rho);*/

    int a = 5;
    
    return 0;
}