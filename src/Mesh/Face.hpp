#ifndef FACE_HPP
#define FACE_HPP

#include <vector>
#include <memory>

#include "Vector3.hpp"

class Face
{
    public:
        enum faceType{GENERAL, TRIANGULAR, QUADRILATERAL};

        Face() : type(GENERAL) {};
        Face(std::vector<int> nodesIdx) : nodesIndex(nodesIdx), type(GENERAL) {};
        Face(std::vector<int> nodesIdx, faceType fType) : nodesIndex(nodesIdx), type(fType) {};

        void update(const std::vector<Vector3>& nodeList);
        bool check() const;
        bool equal(const Face& compFace) const;

        void reverseOrientation();

        virtual ~Face();

        std::vector<int> nodesIndex;

        double area;
        Vector3 normalVector;
        Vector3 midpoint;

    protected:
        const int type;

        Vector3 calculateNormalVector(const std::vector<Vector3>& nodeList);
        Vector3 calculateMidpoint(const std::vector<Vector3>& nodeList) const;

};

#endif // FACE_HPP