#include <cmath>
#include <iostream>
#include <fstream>

#include "SpecialGasThermo.hpp"


/*Vars<3> SpecialGasThermo::updateThermo(const Compressible& data, const Compressible& dataOld) const
{
    //treti clen stara hodnota teploty jako odhad pro nelinearni reseni rovnice
    double rho = data.density();
    double T = tFromRhoE(rho, data.internalEnergy(), dataOld.temperature());

    return Vars<3>({T, p(rho, T), a(rho, T)});
}*/

Vars<3> SpecialGasThermo::updateThermo(const Compressible& data) const
{
    double rho = data.density();
    //double TGuess = ((1.29 - 1.0)*(data.totalEnergy() - 0.5*data.absVelocity2()))/specGasConst; //ideal gas guess
    double T = tFromRhoE(rho, data.internalEnergy(), data.temperature());

    return Vars<3>({T, p(rho, T), a(rho, T)});
}

Compressible SpecialGasThermo::primitiveToConservative(const Vars<5>& primitive) const
{
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

Compressible SpecialGasThermo::stagnationState(double TTot, double pTot) const
{
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return Compressible({rhoTot,
                         0.0,
                         0.0,
                         0.0,
                         rhoTot*(e(rhoTot, TTot))},
                         {TTot, pTot, a(rhoTot, TTot)});
}

/*Compressible SpecialGasThermo::isentropicInlet(double pTot, double TTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
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
}*/

/*Compressible SpecialGasThermo::isentropicInletPressureTemperature(double pTot, double TTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}*/

/*Compressible SpecialGasThermo::isentropicInletPressureDensity(double pTot, double rhoTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    //ideal gas as guess
    double TTot = tFromRhoP(rhoTot, pTot, pTot/(specGasConst*rhoTot));

    return isentropicInlet(pTot, TTot, rhoTot, velocityDirection, stateIn);
}*/

Compressible SpecialGasThermo::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    stateIn.setThermoVar(updateThermo(stateIn));

    double pIn = std::min(stateIn.pressure(), pTot);

    double guessRho = stateIn.density();
    double guessT = stateIn.temperature();

    std::pair<double, double> result = RhoTFromSP(sTot, pIn, guessRho, guessT);

    double rho = result.first;
    double T = result.second;

    double absU2 = std::max(2.0*hTot - 2.0*h(rho, T), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e(rho, T))},
                         {T, pIn, a(rho, T)});
}

std::array<double, 3> SpecialGasThermo::initPressureTemperatureInlet(double pTot, double TTot) const
{
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return std::array<double, 3>({rhoTot, s(rhoTot, TTot), h(rhoTot, TTot)});
}