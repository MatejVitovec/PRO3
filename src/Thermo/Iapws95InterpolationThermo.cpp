#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95InterpolationThermo.hpp"

#include"Interpolation/BiLinearInterpolation.hpp"
//#include"Interpolation/BiQuadraticInterpolation.hpp"


Iapws95InterpolationThermo::Iapws95InterpolationThermo(InterpolationType interpolationType) : Thermo(), Iapws95()
{
    switch (interpolationType)
    {
    case BILINEAR:
        pressureInterpolationFromRhoE = std::make_unique<BiLinearInterpolation>(std::vector<int>{200},
            std::vector<int>{100, 100},
            std::vector<double>{1/1188.87, 1/0.008},
            std::vector<double>{2005000.0, 2650000.0, 4085000.27},
            Interpolation::LOGINV,
            Interpolation::NONE);

        soundSpeedInterpolationFromRhoE = std::make_unique<BiLinearInterpolation>(std::vector<int>{200},
            std::vector<int>{100, 100},
            std::vector<double>{1/1188.87, 1/0.008},
            std::vector<double>{2005000.0, 2650000.0, 4085000.27},
            Interpolation::LOGINV,
            Interpolation::NONE); //TODO

        temperatureInterpolationFromRhoE = std::make_unique<BiLinearInterpolation>(std::vector<int>{200},
            std::vector<int>{100, 100},
            std::vector<double>{1/1188.87, 1/0.008},
            std::vector<double>{2005000.0, 2650000.0, 4085000.27},
            Interpolation::LOGINV,
            Interpolation::NONE); //TODO

        /*energyFromRhoPAux = std::make_unique<BiLinearInterpolation>(std::vector<int>{200},
            std::vector<int>{100, 100},
            std::vector<double>{1/1188.87, 1/0.008},
            std::vector<double>{2005000.0, 2650000.0, 4085000.27},
            Interpolation::LOGINV,
            Interpolation::NONE); //TODO*/

        break;
    case BIQUADRATIC:
        /* code */
        break;
    default:
        std::cout << "Error: Bad interpolation type" << std::endl;
    }

    pressureInterpolationFromRhoE->calcCoeffs([=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
    
    soundSpeedInterpolationFromRhoE->calcCoeffs([=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
    
    temperatureInterpolationFromRhoE->calcCoeffs([=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

    //energyFromRhoPAux->calcCoeffs([=](double density, double energy) { return 0.0; }); //TODO
}

Vars<3> Iapws95InterpolationThermo::updateThermo(const Compressible& data, const Compressible& dataOld) const
{
    double pressure = pressureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());

    return Vars<3>({temperature, pressure, soundSpeed}); //T = 0 - nepotrebuju
}

Vars<3> Iapws95InterpolationThermo::updateThermo(const Compressible& data) const
{
    double pressure = pressureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE->calc(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE->calc(data.density(), data.internalEnergy());

    return Vars<3>({temperature, pressure, soundSpeed}); //T = 0 - nepotrebuju
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

Compressible Iapws95InterpolationThermo::isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const
{
    // FULL IAPWS95
    stateIn.setThermoVar(updateThermo(stateIn, wrOld));

    double pIn = std::min(stateIn.pressure(), pTot);
    double sTot = s(rhoTot, TTot);

    double guessRho = stateIn.density();
    double guessT = wrOld.temperature();

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

Compressible Iapws95InterpolationThermo::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const
{
    // FULL IAPWS95
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn, wrOld);
}

Compressible Iapws95InterpolationThermo::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn, Compressible wrOld) const
{
    // FULL IAPWS95
    double TTot = tFromRhoP(rhoTot, pTot, pTot/(specGasConst*rhoTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn, wrOld);
}
