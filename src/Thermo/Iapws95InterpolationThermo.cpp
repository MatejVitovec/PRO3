#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95InterpolationThermo.hpp"

#include"Interpolation/BiLinearInterpolation.hpp"
#include"Interpolation/BiQuadraticInterpolation.hpp"


Iapws95InterpolationThermo::Iapws95InterpolationThermo(InterpolationType interpolationType) : Thermo(), Iapws95()
{
    std::vector<double> densityGrid = createInterpolationAxis(std::vector<int>{200}, std::vector<double>{1/1188.87, 1/0.008}, Transformation::LOGINV);
    std::vector<double> internalEnergyGrid = createInterpolationAxis(std::vector<int>{100, 100}, std::vector<double>{2005000.0, 2650000.0, 4085000.27}, Transformation::NONE);

    switch (interpolationType)
    {
    case BILINEAR:
        pressureInterpolationFromRhoE    = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
        soundSpeedInterpolationFromRhoE  = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
        temperatureInterpolationFromRhoE = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

        //energyFromRhoPAux = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return 0.0; }); TODO
        break;

    case BIQUADRATIC:
        //TODO derivace - nebo numericky
        pressureInterpolationFromRhoE    = std::make_unique<BiQuadraticInterpolation>(densityGrid, internalEnergyGrid,
            [=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

        soundSpeedInterpolationFromRhoE  = std::make_unique<BiQuadraticInterpolation>(densityGrid, internalEnergyGrid,
            [=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

        temperatureInterpolationFromRhoE = std::make_unique<BiQuadraticInterpolation>(densityGrid, internalEnergyGrid,
            [=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

                std::cout << "BiQuadratic interpolation thermo init ok" << std::endl;

        /*energyFromRhoPAux = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid,
            [=](double density, double energy) { return 0.0; },
            [=](double density, double energy) { return 0.0; },
            [=](double density, double energy) { return 0.0; },
            [=](double density, double energy) { return 0.0; });*/
        break;

    default:
        std::cout << "Error: Bad interpolation type" << std::endl;
    }

}

Vars<3> Iapws95InterpolationThermo::updateThermo(const Compressible& data) const
{
    double pressure = pressureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());

    return Vars<3>({temperature, pressure, soundSpeed});
}

Vars<3> Iapws95InterpolationThermo::updateThermo(const Primitive& data) const
{
    /*double pressure = pressureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());*/

    double pressure = 0.0;
    double soundSpeed = 0.0;
    double temperature = 0.0;

    return Vars<3>({temperature, pressure, soundSpeed});
}

Compressible Iapws95InterpolationThermo::primitiveToConservative(const Vars<5>& primitive) const
{
    /*double density = primitive[0];
    double pressure = primitive[4];

    double energy = pressureInterpolationFromRhoE->calcInverseY(density, pressure, energyFromRhoPAux->calc(density, pressure));
    double soundSpeed = soundSpeedInterpolationFromRhoE->calc(density, energy);

    return Compressible({density,
                         density*primitive[1],
                         density*primitive[2],
                         density*primitive[3],
                         density*(energy + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {0.0, pressure, soundSpeed});*/

    //FULL IAPWS95
    double rho = primitive[0];
    double p = primitive[4];
    double T = tFromRhoP(rho, p, p/(specGasConst*rho)); //odhad T pomoci idealniho plynu

    return Compressible({rho,
                         rho*primitive[1],
                         rho*primitive[2],
                         rho*primitive[3],
                         rho*(e(rho, T) + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {T, p, a(rho, T)});
}

Compressible Iapws95InterpolationThermo::stagnationState(double TTot, double pTot) const
{
    // FULL IAPWS95
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return Compressible({rhoTot,
                         0.0,
                         0.0,
                         0.0,
                         rhoTot*(e(rhoTot, TTot))},
                         {TTot, pTot, a(rhoTot, TTot)});
}

Compressible Iapws95InterpolationThermo::isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    // FULL IAPWS95
    stateIn.setThermoVar(updateThermo(stateIn));

    double pIn = std::min(stateIn.pressure(), pTot);
    double sTot = s(rhoTot, TTot);

    double guessRho = stateIn.density();
    double guessT = stateIn.temperature();

    std::pair<double, double> result = RhoTFromSP(sTot, pIn, guessRho, guessT);

    double rho = result.first;
    double T = result.second;

    double absU2 = std::max(2.0*h(rhoTot, TTot) - 2.0*h(rho, T), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e(rho, T))},
                         {T, pIn, a(rho, T)});
}

Compressible Iapws95InterpolationThermo::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    // FULL IAPWS95
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

Compressible Iapws95InterpolationThermo::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    // FULL IAPWS95
    double TTot = tFromRhoP(rhoTot, pTot, pTot/(specGasConst*rhoTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}

std::vector<double> Iapws95InterpolationThermo::createInterpolationAxis(std::vector<int> gridSize, std::vector<double> boundary, Transformation transformation)
{
    for (int i = 0; i < boundary.size(); i++)
    {
        switch(transformation) 
        {
            case LOG:
                boundary[i] = log(boundary[i]);
                break;
            case LOG10:
                boundary[i] = log10(boundary[i]);
                break;
            case LOGINV:
                boundary[i] = log(1/boundary[i]);
                break;
        }
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

    auto backTransform = [transformation](double x)
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
    };

    int idx = 0;
    out[idx] = backTransform(boundary[0]);
    idx++;

    for (int ii = 0; ii < gridSize.size(); ii++)
    {
        for (int i = 0; i < gridSize[ii]; i++)
        {
            out[idx] = backTransform(boundary[ii] + dx[ii]*(i+1));
            idx++;
        }        
    }

    if(out[0] > out[1])
    {
        std::reverse(out.begin(), out.end());
    }

    return out;
}