#ifndef IAPWS95INTERPOLATIONTHERMO
#define IAPWS95INTERPOLATIONTHERMO

#include <memory>

#include "Thermo.hpp"
#include "StateEquations/Iapws95.hpp"
#include "Interpolation/Interpolation.hpp"
#include "StateEquations/NonLinearSolver/NewtonMethod.hpp"

template <typename INTERPOLATION>
class Iapws95InterpolationThermo : public Thermo, Iapws95
{
    public:
        enum InterpolationType{BILINEAR, BIQUADRATIC};

        Iapws95InterpolationThermo();

        Vars<3> updateThermo(const Compressible& data) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;


    private:
        INTERPOLATION pressureInterpolationFromRhoE;
        INTERPOLATION soundSpeedInterpolationFromRhoE;
        INTERPOLATION temperatureInterpolationFromRhoE;
        INTERPOLATION entropyInterpolationFromRhoE;

        INTERPOLATION energyInterpolationFromRhoPAux;

        NewtonMethod newtonSolver;
};

template <typename INTERPOLATION>
Iapws95InterpolationThermo<INTERPOLATION>::Iapws95InterpolationThermo() : Thermo(), Iapws95(), newtonSolver(NewtonMethod())
{
    std::vector<int> densityGridSize({200});
    std::vector<double> densityBoundary({1/1188.87, 1/0.008});
    Interpolation::Transformation densityTransformation(Interpolation::LOGINV);

    std::vector<int> internalEnergyGridSize({50, 100, 100});
    std::vector<double> internalEnergyBoundary({1800000.0, 2005000.0, 2650000.0, 4085000.27});
    Interpolation::Transformation internalEnergyTransformation(Interpolation::NONE);

    std::vector<int> pressureGridSize({800});
    std::vector<double> pressureBoundary({7000.0, 200000.0});
    Interpolation::Transformation pressureTransformation(Interpolation::LOGINV);

    pressureInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return p(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

    soundSpeedInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return a(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

    temperatureInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

    energyInterpolationFromRhoPAux = INTERPOLATION(
        densityGridSize, pressureGridSize, densityBoundary, pressureBoundary, densityTransformation, pressureTransformation,
        [=](double density, double pressure) { return e(density, tFromRhoP(density, pressure, pressure/(density*461.51805))); });

    entropyInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return s(density, tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
}

template <typename INTERPOLATION>
Vars<3> Iapws95InterpolationThermo<INTERPOLATION>::updateThermo(const Compressible& data) const
{
    double pressure = pressureInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());

    return Vars<3>({temperature, pressure, soundSpeed});
}


template <typename INTERPOLATION>
Compressible Iapws95InterpolationThermo<INTERPOLATION>::primitiveToConservative(const Vars<5>& primitive) const
{
    double density = primitive[0];
    double pressure = primitive[4];

    double energyAux = energyInterpolationFromRhoPAux.calcFastFind(density, pressure);
    
    double energy = pressureInterpolationFromRhoE.calcInverseY(density, pressure, energyAux);
    double soundSpeed = soundSpeedInterpolationFromRhoE.calcFastFind(density, energy);

    return Compressible({density,
                         density*primitive[1],
                         density*primitive[2],
                         density*primitive[3],
                         density*(energy + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))},
                         {0.0, pressure, soundSpeed});
}

template <typename INTERPOLATION>
Compressible Iapws95InterpolationThermo<INTERPOLATION>::stagnationState(double TTot, double pTot) const
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

template <typename INTERPOLATION>
Compressible Iapws95InterpolationThermo<INTERPOLATION>::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn) const
{
    double pIn = std::min(stateIn.pressure(), pTot);

    double guessRho = stateIn.density();
    double guess_e = stateIn.internalEnergy();

    std::pair<double, double> result = newtonSolver.solve([=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFind(val1, val2)/sTot - 1; },
                                                          [=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFindDiffX(val1, val2)/sTot; },
                                                          [=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFindDiffY(val1, val2)/sTot; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFind(val1, val2)/pIn - 1; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFindDiffX(val1, val2)/pIn; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFindDiffY(val1, val2)/pIn; },
                                                          guessRho,
                                                          guess_e);

    double rho = result.first;
    double e = result.second;

    double absU2 = std::max(2.0*hTot - 2.0*(e + pressureInterpolationFromRhoE.calcFastFind(rho, e)/rho), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e)},
                        {temperatureInterpolationFromRhoE.calcFastFind(rho, e),
                         pIn,
                         soundSpeedInterpolationFromRhoE.calcFastFind(rho, e)});
}

template <typename INTERPOLATION>
std::array<double, 3> Iapws95InterpolationThermo<INTERPOLATION>::initPressureTemperatureInlet(double pTot, double TTot) const
{
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return std::array<double, 3>({rhoTot, s(rhoTot, TTot), h(rhoTot, TTot)});
}

#endif // IAPWS95INTERPOLATIONTHERMO