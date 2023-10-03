#include <cmath>
#include <iostream>
#include "Cell.hpp"
#include "Vector3.hpp"


void Cell::update(const std::vector<Face>& faceList)
{
    double vol = 0.0;
    Vector3 cen = Vector3(0.0, 0.0, 0.0);
    Vector3 projArea = Vector3(0.0, 0.0, 0.0);

    for (auto & faceIndex : ownFaceIndex)
    {
        Face face = faceList[faceIndex];
        Vector3 midpoint = face.midpoint;
        Vector3 normal = face.area*face.normalVector;
        vol += dot(midpoint, normal);
        cen = cen + std::pow(norm2(midpoint), 2.0)*normal;
        projArea = projArea + abs(normal);
    }

    for (auto & faceIndex : neighborFaceIndex)
    {
        Face face = faceList[faceIndex];
        Vector3 midpoint = face.midpoint;
        Vector3 normal = (-1.0)*face.area*face.normalVector;
        vol += dot(midpoint, normal);
        cen = cen + std::pow(norm2(midpoint), 2.0)*normal;
        projArea = projArea + abs(normal);
    }
    
    volume = vol/3.0;
    center = (0.5/volume)*cen;
    projectedArea = 0.5*projArea;
}

int Cell::getVtkType() const
{
    if(type == TETRAHEDRON) return 10;
    if(type == HEXAHEDRON) return 12;
    if(type == PRISM) return 13;
    if(type == PYRAMID) return 14;
    return 0;
}

std::vector<Face> Cell::createFaces()
{
    switch(type)
    {
        case TETRAHEDRON:
            return std::vector<Face>{Face(std::vector<int>{nodesIndex[0], nodesIndex[2], nodesIndex[1]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[1], nodesIndex[3]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[3], nodesIndex[2]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[1], nodesIndex[2], nodesIndex[3]}, Face::TRIANGULAR)};

        case HEXAHEDRON:
            return std::vector<Face>{Face(std::vector<int>{nodesIndex[0], nodesIndex[1], nodesIndex[5], nodesIndex[4]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[2], nodesIndex[3], nodesIndex[7], nodesIndex[6]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[4], nodesIndex[7], nodesIndex[3]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[4], nodesIndex[5], nodesIndex[6], nodesIndex[7]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[1], nodesIndex[2], nodesIndex[6], nodesIndex[5]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[3], nodesIndex[2], nodesIndex[1]}, Face::QUADRILATERAL)};

        case PRISM:
            return std::vector<Face>{Face(std::vector<int>{nodesIndex[0], nodesIndex[2], nodesIndex[1]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[3], nodesIndex[4], nodesIndex[5]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[1], nodesIndex[4], nodesIndex[3]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[3], nodesIndex[5], nodesIndex[2]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[1], nodesIndex[2], nodesIndex[5], nodesIndex[4]}, Face::QUADRILATERAL)};


        case PYRAMID:
            return std::vector<Face>{Face(std::vector<int>{nodesIndex[0], nodesIndex[3], nodesIndex[2], nodesIndex[1]}, Face::QUADRILATERAL),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[1], nodesIndex[4]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[1], nodesIndex[2], nodesIndex[4]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[2], nodesIndex[3], nodesIndex[4]}, Face::TRIANGULAR),
                                     Face(std::vector<int>{nodesIndex[0], nodesIndex[4], nodesIndex[3]}, Face::TRIANGULAR)};
        
        default:
            std::cout << "ERROR - Nedefinovaný typ buňky" << std::endl;
            return std::vector<Face>();
    }
}

Cell::~Cell()
{
    
}