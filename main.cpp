#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>

#include "src/Thermo/Iapws95Thermo.hpp"



int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    std::shared_ptr<Iapws95Thermo> thermo = std::make_shared<Iapws95Thermo>();

    double M2is = thermo->calcM2is(435000.0, 350.0, 163000.0);

    std::cout << "M2is: " << M2is << std::endl;

    return 0;
}