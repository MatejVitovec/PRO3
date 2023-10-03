#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <memory>
#include <string>

#include "Vector3.hpp"
#include "Cell.hpp"
#include "Face.hpp"
#include "Boundary.hpp"


class Mesh
{
    public:
        Mesh() {};

        const std::vector<Vector3>& getNodeList() const;
        const std::vector<Cell>& getCellList() const;
        const std::vector<Face>& getFaceList() const;
        const std::vector<Boundary>& getBoundaryList() const;
        const std::vector<int>& getOwnerIndexList() const;
        const std::vector<int>& getNeighborIndexList() const;

        int getNodesSize() const;
        int getFacesSize() const;
        int getCellsSize() const;

        void update();
        void loadGmsh2(std::string fileName);

    private:
        void createFaces();
        void createBoundaryFacesGmsh(const std::vector<std::vector<std::string>>& physicalNamesGmsh, const std::vector<std::vector<std::string>>& elementsGmsh);
        
        void updateCells();
        void updateFaces();

        bool checkFaces() const;

        std::vector<Vector3> nodeList;
        std::vector<Cell> cellList;
        std::vector<Face> faceList;
        std::vector<Boundary> boundaryList;
        
        std::vector<int> ownerIndexList;
        std::vector<int> neighborIndexList;

        //load GMSH
        std::vector<std::string> readFile(std::string fileName);
        std::vector<std::vector<std::string>> parseBlockDataGmsh(const std::vector<std::string>& dataIn, std::string blockName);
        void createNodesGmsh(const std::vector<std::vector<std::string>>& nodesGmsh);
        void createCellsGmsh(const std::vector<std::vector<std::string>>& elementsGmsh);
        void createBoundariesGmsh(const std::vector<std::vector<std::string>>& physicalNamesGmsh, const std::vector<std::vector<std::string>>& elementsGmsh);

};

#endif // MESH_HPP