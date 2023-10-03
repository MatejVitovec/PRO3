#ifndef OUTPUTCFD_HPP
#define OUTPUTCFD_HPP

#include "Mesh/Mesh.hpp"
#include "Field.hpp"
#include "Compressible.hpp"
#include <string>

void outputVTK(std::string filename, const Mesh& m, const Field<Compressible>& u);

void saveResidual(std::string filename, Vars<5> res);

#endif //OUTPUTCFD_HPP