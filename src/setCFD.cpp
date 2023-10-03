#include "setCFD.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "BoundaryCondition/PressureTemperatureInlet.hpp"
#include "BoundaryCondition/PressureDensityInlet.hpp"
#include "BoundaryCondition/PressureOutlet.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"
#include "BoundaryCondition/Wall.hpp"



std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryCondition(const Mesh& mesh)
{
	std::vector<std::unique_ptr<BoundaryCondition>> out;

	std::vector<std::string> boundaryNames;
	std::vector<int> boundaryTypes;

	boundaryNames.push_back("inlet");
	boundaryNames.push_back("outlet");
	boundaryNames.push_back("top");
	boundaryNames.push_back("bottom");
	boundaryNames.push_back("front");
	boundaryNames.push_back("back");

	boundaryTypes.push_back(BoundaryCondition::PRESSUREDENSITYINLET);
	boundaryTypes.push_back(BoundaryCondition::PRESSUREOUTLET);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);

	const std::vector<Boundary>& meshBoundaryList = mesh.getBoundaryList();
	
	for (int i = 0; i < boundaryTypes.size(); i++)
	{
		Boundary aux;
		bool exist = false;
		for (const auto & meshBoundary : meshBoundaryList)
		{
			if(meshBoundary.boundaryConditionName == boundaryNames[i])
			{
				aux = meshBoundary;
				exist = true;
				break;
			}
		}

		if(!exist)
		{
			std::cout << "Chyba BC" << std::endl;
			continue;
		}

		if(boundaryTypes[i] == BoundaryCondition::PRESSURETEMPERATUREINLET)
		{
			out.push_back(std::make_unique<PressureTemperatureInlet>(aux, 0.7143, 0.5102, Vars<3>({1.0, 0.0, 0.0})));
		}
		else if(boundaryTypes[i] == BoundaryCondition::PRESSUREDENSITYINLET)
		{
			out.push_back(std::make_unique<PressureDensityInlet>(aux, 1.0/1.4, 1.0, Vars<3>({1.0, 0.0, 0.0})));
		}
		else if(boundaryTypes[i] == BoundaryCondition::PRESSUREOUTLET)
		{
			out.push_back(std::make_unique<PressureOutlet>(aux, 0.5264));
		}
		else if(boundaryTypes[i] == BoundaryCondition::FREEBOUNDARY)
		{
			out.push_back(std::make_unique<FreeBoundary>(aux));
		}
		else if(boundaryTypes[i] == BoundaryCondition::WALL)
		{
			out.push_back(std::make_unique<Wall>(aux));
		}
	}

	return out;

}