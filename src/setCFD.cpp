#include "setCFD.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "BoundaryCondition/PressureTemperatureInlet.hpp"
#include "BoundaryCondition/PressureDensityInlet.hpp"
#include "BoundaryCondition/PressureOutlet.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"
#include "BoundaryCondition/Wall.hpp"
#include "BoundaryCondition/Periodicity.hpp"



std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryCondition(const Mesh& mesh)
{
	std::vector<std::unique_ptr<BoundaryCondition>> out;

	std::vector<std::string> boundaryNames;
	std::vector<int> boundaryTypes;

	/*boundaryNames.push_back("inlet");
	boundaryNames.push_back("outlet");
	boundaryNames.push_back("top");
	boundaryNames.push_back("bottom");
	boundaryNames.push_back("front");
	boundaryNames.push_back("back");
	boundaryTypes.push_back(BoundaryCondition::PRESSURETEMPERATUREINLET);
	boundaryTypes.push_back(BoundaryCondition::PRESSUREOUTLET);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);*/

	boundaryNames.push_back("inlet");
	boundaryNames.push_back("outlet");
	boundaryNames.push_back("wall");
	boundaryNames.push_back("wall2");
	boundaryNames.push_back("periodbeg");
	boundaryNames.push_back("periodend");
	

	boundaryTypes.push_back(BoundaryCondition::PRESSURETEMPERATUREINLET);
	boundaryTypes.push_back(BoundaryCondition::PRESSUREOUTLET);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::WALL);
	boundaryTypes.push_back(BoundaryCondition::PERIODICITY);
	boundaryTypes.push_back(100);

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
			//out.push_back(std::make_unique<PressureTemperatureInlet>(aux, 99980.0, 422.65, Vars<3>({1.0, 0.0, 0.0})));
			out.push_back(std::make_unique<PressureTemperatureInlet>(aux, 98071.7, 298.65, vector3toVars(angleAngleToUnit(19.3, 0.0))));
		}
		else if(boundaryTypes[i] == BoundaryCondition::PRESSUREDENSITYINLET)
		{
			out.push_back(std::make_unique<PressureDensityInlet>(aux, 1.0/1.4, 1.0, Vars<3>({1.0, 0.0, 0.0})));
		}
		else if(boundaryTypes[i] == BoundaryCondition::PRESSUREOUTLET)
		{
			//out.push_back(std::make_unique<PressureOutlet>(aux, /*0.5264*/80000.0));
			out.push_back(std::make_unique<PressureOutlet>(aux, 40548.109));
		}
		else if(boundaryTypes[i] == BoundaryCondition::FREEBOUNDARY)
		{
			out.push_back(std::make_unique<FreeBoundary>(aux));
		}
		else if(boundaryTypes[i] == BoundaryCondition::WALL)
		{
			out.push_back(std::make_unique<Wall>(aux));
		}
		else if(boundaryTypes[i] == BoundaryCondition::PERIODICITY)
		{
			out.push_back(std::make_unique<Periodicity>(aux, Vector3(0.0, 0.0551168, 0.0), "aaa", mesh));
		}
		else if(boundaryTypes[i] == 100)
		{
			out.push_back(std::make_unique<Periodicity>(aux, Vector3(0.0, -0.0551168, 0.0), "bbb", mesh));
		}
	}

	return out;

}