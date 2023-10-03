#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <vector>
#include <memory>
#include <string>

class Boundary
{
    public:
        Boundary() : boundaryConditionName("empty") {};
        Boundary(std::string name) : boundaryConditionName(name) {};
        //Boundary(std::vector<int> fIndexes) : facesIndex(fIndexes) {};

        //std::vector<int> nodesIndex;
        std::vector<int> facesIndex;

        std::string boundaryConditionName;

    private:

};

#endif // BOUNDARY_HPP