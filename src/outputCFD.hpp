#ifndef OUTPUTCFD_HPP
#define OUTPUTCFD_HPP

#include "Mesh/Mesh.hpp"
#include "Field.hpp"
#include "Compressible.hpp"
#include "Mat.hpp"
#include <string>

namespace outputCFD
{
    void outputVTK(std::string fileName, const Mesh& mesh, const Field<Compressible>& w);

    void outputVTKPeriodicBoundary(std::string fileName, const Mesh& m, const Field<Compressible>& w, Vector3 shift);

    void saveData(std::string fileName, const Field<Compressible>& w);

    void saveResidual(std::string fileName, Vars<5> res);
    void saveValue(std::string fileName, double val);

    void saveFieldOnBoundary(std::string fileName, std::string boundaryName, const Mesh& mesh, const Field<Compressible>& w);

    void saveLimiters(Field<Vars<5>> phi, const Mesh& mesh);
    void saveGradients(Field<Mat<5,3>> grad, const Mesh& mesh);
}

#endif //OUTPUTCFD_HPP