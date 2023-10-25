#ifndef SETCFD_HPP
#define SETCFD_HPP

#include <string>

#include "Mesh/Mesh.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"


std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryCondition(const Mesh& mesh);


#endif //SETCFD_HPP