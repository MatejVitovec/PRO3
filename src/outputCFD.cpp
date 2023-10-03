#include "outputCFD.hpp"
#include <fstream>
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

void outputVTK(std::string filename, const Mesh& mesh, const Field<Compressible>& w)
{
    const std::vector<Vector3>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();

    int cellSize = mesh.getCellsSize();

	std::ofstream f;
	f.open(filename, std::ios::out);
	
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

 	//f << "SCALARS |u| float\n"; 
	//f << "LOOKUP_TABLE default\n";

    //for (int i=0; i<m.nc; i++) {
	//	f << u[i].u().norm() << "\n";
	//}		

	//f << "VECTORS u float\n"; 
 	f << "SCALARS u float 3\n"; 
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
	
 	f << "SCALARS x float 3\n"; 
	f << "LOOKUP_TABLE default\n";

    for (int i = 0; i < cellSize; i++)
    {
		f << cellList[i].center << "\n";
	}
	f << std::endl;

	f.close();
}


void saveResidual(std::string filename, Vars<5> res)
{

	std::ofstream f;
	f.open(filename, std::ios_base::app);

	//f << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << std::endl;
	f << res[0] <<std::endl;

	f.close();
}