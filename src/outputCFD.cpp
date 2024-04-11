#include "outputCFD.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

int calculateCellNodeSize(const Mesh& mesh)
{
    const std::vector<Cell>& cellList = mesh.getCellList();
    int num = 0;

    for (int i = 0; i < cellList.size(); i++)
    {
        num += (1 + cellList[i].nodesIndex.size());
    }
    
    return num;
}

double roundToZero(double in)
{
	if(in < 10e-12 && in > -10e-12)
	{
		return 0.0;
	}
	return in;
}

void outputCFD::outputVTK(std::string fileName, const Mesh& mesh, const Field<Compressible>& w)
{
    const std::vector<Vector3>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();

    int cellSize = mesh.getCellsSize();

	std::ofstream f;
	f.open(fileName, std::ios::out);
	
	f << "# vtk DataFile Version 1.0\n";
	f << "unstructured grid\n";
	f << "ascii\n";
	f << "DATASET UNSTRUCTURED_GRID\n";
	
	f << "points " << mesh.getNodesSize() << " float\n";
	
	for (int i = 0; i < mesh.getNodesSize(); i++)
    {
		f << nodeList[i] << "\n";
	}
	
	f << "cells " << cellSize << " " << calculateCellNodeSize(mesh) << "\n";
	
	for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i] << "\n";
	}
	
	f << "cell_types " << cellSize << "\n";
	
	for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i].getVtkType() << "\n";
	}
	
	f << "CELL_DATA " << cellSize << "\n";
 	f << "SCALARS rho float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].density()) << "\n";
	}

	//f << "VECTORS u float\n"; 
 	f << "SCALARS U float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].velocityU()) << " " << roundToZero(w[i].velocityV()) << " " << roundToZero(w[i].velocityW()) << "\n";
	}	

 	f << "SCALARS e float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].internalEnergy()) << "\n";
	}
	
	f << "SCALARS p float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].pressure()) << "\n";
	}

	f << "SCALARS M float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].machNumber()) << "\n";
	}

	f << "SCALARS T float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].temperature()) << "\n";
	}
	
	f << std::endl;

	f.close();
}

void outputCFD::saveData(std::string fileName, const Field<Compressible>& w)
{
	std::ofstream f;
	f.open(fileName, std::ios::out);

	for (int i = 0; i < w.size(); i++)
	{
		for (int j = 0; j < 5; j++)
		{
			f << roundToZero(w[i][j]) << " ";
		}
		f << "\n";
	}

	f << std::endl;

	f.close();
}


void outputCFD::saveResidual(std::string fileName, Vars<5> res)
{

	std::ofstream f;
	f.open(fileName, std::ios_base::app);

	//f << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << std::endl;
	f << res[0] <<std::endl;

	f.close();
}

void outputCFD::saveValue(std::string fileName, double val)
{
	std::ofstream f;
	f.open(fileName, std::ios_base::app);
	f << val << std::endl;
	f.close();
}

void outputCFD::saveFieldOnBoundary(std::string fileName, std::string boundaryName, const Mesh& mesh, const Field<Compressible>& w)
{
	const std::vector<Boundary>& boundary = mesh.getBoundaryList();
	const std::vector<Face>& faces = mesh.getFaceList();
	const std::vector<int>& owners = mesh.getOwnerIndexList();

	int boundaryIndex = 0;
	bool found = false;

	for (boundaryIndex = 0; boundaryIndex < boundary.size(); boundaryIndex++)
	{
		if (boundary[boundaryIndex].boundaryConditionName == boundaryName)
		{
			found = true;
			break;
		}		
	}
	if(!found)
	{
		std::cout << "Error: spatne zadany nazev BC pro vypsani hodnot" << std::endl;
	}
	
	std::ofstream f;
	f.open(fileName, std::ios::out);
	//f << "#x y z p\n";

	for (int i = 0; i < boundary[boundaryIndex].facesIndex.size(); i++)
	{
		int faceIndex = boundary[boundaryIndex].facesIndex[i];

		//f << faces[faceIndex].midpoint.x << " " << faces[faceIndex].midpoint.y << " " << faces[faceIndex].midpoint.z << " " << w[owners[faceIndex]].pressure() << "\n";
		f << faces[faceIndex].midpoint.x << " " << w[owners[faceIndex]].pressure() << "\n";
	}
	
}


/////////////////////////

void outputCFD::outputVTKPeriodicBoundary(std::string fileName, const Mesh& mesh, const Field<Compressible>& w, Vector3 shift)
{
    const std::vector<Vector3>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();

    int cellSize = mesh.getCellsSize();

	std::vector<Cell> shiftedCells = cellList;
	int nodesListSize = nodeList.size();
	for (int i = 0; i < cellList.size(); i++)
	{
		for (int j = 0; j < cellList[i].nodesIndex.size(); j++)
		{
			shiftedCells[i].nodesIndex[j] = cellList[i].nodesIndex[j] + nodesListSize;
		}
	}


	/*const std::vector<Face>& faceList = mesh.getFaceList();
    double minFaceSize = 100000.0;
    for (int i = 0; i < faceList.size(); i++)
    {
        if (minFaceSize > faceList[i].area)
        {
            minFaceSize = faceList[i].area;
        }        
    }
    double numTol = minFaceSize/2.0;*/
	double numTol = 0.000001;


	std::vector<int> duplicityNodes;
	std::vector<int> duplicityNodesNew;

	for (int i = 0; i < nodeList.size(); i++)
	{
		Vector3  node = nodeList[i] - shift;

		for (int j = 0; j < nodeList.size(); j++)
		{
			if (norm2(node - nodeList[j]) < numTol)
			{
				duplicityNodes.push_back(j);
				duplicityNodesNew.push_back(i);
				break;
			}
		}
	}
	
	for (int i = duplicityNodes.size() - 1; i >= 0; i--)
	{
		for (int j = 0; j < shiftedCells.size(); j++)
		{
			for (int k = 0; k < shiftedCells[j].nodesIndex.size(); k++)
			{
				if ((duplicityNodes[i] + nodesListSize) == shiftedCells[j].nodesIndex[k])
				{
					shiftedCells[j].nodesIndex[k] = duplicityNodesNew[i];
					duplicityNodes.pop_back();
					duplicityNodesNew.pop_back();
				}
			}
		}
	}

	std::ofstream f;
	f.open(fileName, std::ios::out);
	
	f << "# vtk DataFile Version 1.0\n";
	f << "unstructured grid\n";
	f << "ascii\n";
	f << "DATASET UNSTRUCTURED_GRID\n";
	
	f << "points " << 2*mesh.getNodesSize() << " float\n";
	
	for (int i = 0; i < mesh.getNodesSize(); i++)
    {
		f << nodeList[i] << "\n";
	}
	for (int i = 0; i < mesh.getNodesSize(); i++)
    {
		f << (nodeList[i] + shift) << "\n";
	}
	
	f << "cells " << 2*cellSize << " " << 2*calculateCellNodeSize(mesh) << "\n";
	
	for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i] << "\n";
	}	

	for (int i = 0; i < cellSize; i++)
    {
		f << shiftedCells[i] << "\n";
	}
	
	f << "cell_types " << 2*cellSize << "\n";
	
	for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i].getVtkType() << "\n";
	}

	for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i].getVtkType() << "\n";
	}
	
	f << "CELL_DATA " << 2*cellSize << "\n";
 	f << "SCALARS rho float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].density()) << "\n";
	}
	for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].density()) << "\n";
	}

 	f << "SCALARS U float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].velocityU()) << " " << roundToZero(w[i].velocityV()) << " " << roundToZero(w[i].velocityW()) << "\n";
	}
	for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].velocityU()) << " " << roundToZero(w[i].velocityV()) << " " << roundToZero(w[i].velocityW()) << "\n";
	}	

 	f << "SCALARS e float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].internalEnergy()) << "\n";
	}
	for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].internalEnergy()) << "\n";
	}
	
	f << "SCALARS p float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].pressure()) << "\n";
	}
	for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].pressure()) << "\n";
	}

	f << "SCALARS M float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].machNumber()) << "\n";
	}
	for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].machNumber()) << "\n";
	}

	f << "SCALARS T float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].temperature()) << "\n";
	}
    for (int i = 0; i < cellSize; i++)
    {
		f << roundToZero(w[i].temperature()) << "\n";
	}
	
	f << std::endl;

	f.close();
}

void outputCFD::saveLimiters(Field<Vars<5>> phi, const Mesh& mesh)
{
	const std::vector<Cell>& cellList = mesh.getCellList();

	std::ofstream f;
	f.open("../results/limiters.txt", std::ios::out);

	f << "SCALARS phiRho float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(phi[i][0]) << "\n";
	}

	f << "SCALARS phiU float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(phi[i][1]) << "\n";
	}

	f << "SCALARS phiV float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(phi[i][2]) << "\n";
	}

	f << "SCALARS phiW float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(phi[i][3]) << "\n";
	}

	f << "SCALARS phiE float\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(phi[i][4]) << "\n";
	}


	f << std::endl;

	f.close();
}

void outputCFD::saveGradients(Field<Mat<5,3>> grad, const Mesh& mesh)
{
	const std::vector<Cell>& cellList = mesh.getCellList();

	std::ofstream f;
	f.open("../results/rhoGrad.txt", std::ios::out);

 	f << "SCALARS rhoGrad float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(grad[i][0][0]) << " " << roundToZero(grad[i][0][1]) << " " << roundToZero(grad[i][0][2]) << "\n";
	}

	f << std::endl;

	f.close();

	f.open("../results/rhoUGrad.txt", std::ios::out);

 	f << "SCALARS rhoUGrad float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(grad[i][1][0]) << " " << roundToZero(grad[i][1][1]) << " " << roundToZero(grad[i][1][2]) << "\n";
	}

	f << std::endl;

	f.close();

	f.open("../results/rhoVGrad.txt", std::ios::out);

 	f << "SCALARS rhoVGrad float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(grad[i][2][0]) << " " << roundToZero(grad[i][2][1]) << " " << roundToZero(grad[i][2][2]) << "\n";
	}

	f << std::endl;

	f.close();

	f.open("../results/rhoWGrad.txt", std::ios::out);

 	f << "SCALARS rhoWGrad float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(grad[i][3][0]) << " " << roundToZero(grad[i][3][1]) << " " << roundToZero(grad[i][3][2]) << "\n";
	}

	f << std::endl;

	f.close();

	f.open("../results/rhoEGrad.txt", std::ios::out);

 	f << "SCALARS rhoEGrad float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellList.size(); i++)
    {
		f << roundToZero(grad[i][4][0]) << " " << roundToZero(grad[i][4][1]) << " " << roundToZero(grad[i][4][2]) << "\n";
	}

	f << std::endl;

	f.close();
}