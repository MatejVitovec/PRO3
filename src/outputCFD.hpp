#ifndef OUTPUTCFD_HPP
#define OUTPUTCFD_HPP

#include "Mesh/Mesh.hpp"
#include "Field.hpp"
#include "Compressible.hpp"
#include <string>

namespace outputCFD
{
    void outputVTK(std::string fileName, const Mesh& m, const Field<Compressible>& u);

    void saveResidual(std::string fileName, Vars<5> res);

    void saveFieldOnBoundary(std::string fileName, std::string boundaryName, const Mesh& mesh, const Field<Compressible>& w);
}

#endif //OUTPUTCFD_HPP