#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Iapws95InterpolationThermo.hpp"

#include"Interpolation/BiLinearInterpolation.hpp"
#include"Interpolation/BiQuadraticInterpolation.hpp"


Iapws95InterpolationThermo::Iapws95InterpolationThermo(InterpolationType interpolationType) : Thermo(), Iapws95()
{
    std::vector<int> densityGridSize({200});
    std::vector<double> densityBoundary({1/1188.87, 1/0.008});
    Interpolation::Transformation densityTransformation(Interpolation::LOGINV);

    std::vector<int> internalEnergyGridSize({100, 100});
    std::vector<double> internalEnergyBoundary({2005000.0, 2650000.0, 4085000.27});
    Interpolation::Transformation internalEnergyTransformation(Interpolation::NONE);

    std::vector<int> pressureGridSize({800});
    std::vector<double> pressureBoundary({7000.0, 150000.0});
    Interpolation::Transformation pressureTransformation(Interpolation::LOGINV);

    switch (interpolationType)
    {
    case BILINEAR:
        //pressureInterpolationFromRhoE    = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
        //soundSpeedInterpolationFromRhoE  = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
        //temperatureInterpolationFromRhoE = std::make_unique<BiLinearInterpolation>(densityGrid, internalEnergyGrid, [=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

        //energyInterpolationFromRhoPAux = std::make_unique<BiLinearInterpolation>(densityGrid, pressureGridAux, [=](double density, double pressure) { return e(density, tFromRhoP(density, pressure, pressure/(density*461.51805))); });
        
        break;

    case BIQUADRATIC:
        pressureInterpolationFromRhoE    = std::make_unique<BiQuadraticInterpolation>(
            densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
            [=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

        soundSpeedInterpolationFromRhoE  = std::make_unique<BiQuadraticInterpolation>(
            densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
            [=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

        temperatureInterpolationFromRhoE = std::make_unique<BiQuadraticInterpolation>(
            densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
            [=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

        energyInterpolationFromRhoPAux   = std::make_unique<BiQuadraticInterpolation>(
            densityGridSize, pressureGridSize, densityBoundary, pressureBoundary, densityTransformation, pressureTransformation,
            [=](double density, double pressure) { return e(density, tFromRhoP(density, pressure, pressure/(density*461.51805))); });

        std::cout << "BiQuadratic interpolation thermo init ok" << std::endl;

        break;

    default:
        std::cout << "Error: Bad interpolation type" << std::endl;
    }

}

Vars<3> Iapws95InterpolationThermo::updateThermo(const Compressible& data) const
{
    double pressure = pressureInterpolationFromRhoE->calcFastFind(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE->calcFastFind(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE->calcFastFind(data.density(), data.internalEnergy());

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
    double density = primitive[0];
    double pressure = primitive[4];

    double energyAux = energyInterpolationFromRhoPAux->calcFastFind(density, pressure);

    //double Taux = tFromRhoP(density, pressure, pressure/(specGasConst*density)); //odhad T pomoci idealniho plynu
    //double energyAux2 = e(density, Taux);
    //std::cout << "e interpol: " << energyAux << ", e IAPWS: " << energyAux2 << ", p: " << pressure << "rho : " << density << std::endl;
    
    double energy = pressureInterpolationFromRhoE->calcInverseY(density, pressure, energyAux);
    double soundSpeed = soundSpeedInterpolationFromRhoE->calcFastFind(density, energy);

    return Compressible({density,
                         density*primitive[1],
                         density*primitive[2],
                         density*primitive[3],
                         density*(energy + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {0.0, pressure, soundSpeed});

    //FULL IAPWS95
    /*double rho = primitive[0];
    double p = primitive[4];
    double T = tFromRhoP(rho, p, p/(specGasConst*rho)); //odhad T pomoci idealniho plynu

    return Compressible({rho,
                         rho*primitive[1],
                         rho*primitive[2],
                         rho*primitive[3],
                         rho*(e(rho, T) + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {T, p, a(rho, T)});*/
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