#include <cmath>
#include <algorithm>
#include "Face.hpp"


void Face::update(const std::vector<Vector3>& nodeList)
{
    midpoint = calculateMidpoint(nodeList);

    Vector3 normalVectorScale = calculateNormalVector(nodeList);

    area = norm2(normalVectorScale);
    normalVector = normalVectorScale/area;
}

bool Face::check() const
{
    if(std::adjacent_find(nodesIndex.begin(), nodesIndex.end()) != nodesIndex.end())
    {
        return true;
    }

    return false;
}

bool Face::equal(const Face& compFace) const
{
    bool isEqual = true;

    for (auto & compNodeIndex : compFace.nodesIndex)
    {
        if(std::find(nodesIndex.begin(), nodesIndex.end(), compNodeIndex) == nodesIndex.end())
        {
            isEqual = false;
        }
    }
    
    return isEqual;
}

void Face::reverseOrientation()
{
    std::reverse(nodesIndex.begin(), nodesIndex.end());
    //update();
}

Vector3 Face::calculateNormalVector(const std::vector<Vector3>& nodeList)
{
    Vector3 surface = Vector3();
    //Vector3 center = Vector3();

    int i;
    for (i = 0; i < nodesIndex.size() - 1; i++)
    {
        Vector3 auxSurface = cross(nodeList[nodesIndex[i+1]] - nodeList[nodesIndex[i]], (midpoint - nodeList[nodesIndex[i]])/2.0);
        surface = surface + auxSurface;
        //center = center + norm2(auxSurface)*((nodeList[nodesIndex[i+1]] + nodeList[nodesIndex[i]] + midpoint)/3.0);
    }
    //i++;
    Vector3 auxSurface = cross(nodeList[nodesIndex[0]] - nodeList[nodesIndex[i]], (midpoint - nodeList[nodesIndex[i]])/2.0);
    surface = surface + auxSurface;
    //center = center + norm2((nodeList[nodesIndex[0]] + nodeList[nodesIndex[i]] + midpoint)/3.0)*auxSurface;

    //neni midpoint predelat na faceCenter
    //midpoint = center/norm2(surface);

    return surface;
}

Vector3 Face::calculateMidpoint(const std::vector<Vector3>& nodeList) const
{
    Vector3 aux = Vector3();

    int i;
    for (i = 0; i < nodesIndex.size(); i++)
    {
        aux = aux + nodeList[nodesIndex[i]];
    }
    
    return aux / ((double) i);
}

Face::~Face()
{
    
}