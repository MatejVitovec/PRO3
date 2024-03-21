#ifndef GRADIENTSCHEME_HPP
#define GRADIENTSCHEME_HPP

#include "../Mesh/Mesh.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"
#include "../BoundaryCondition/Periodicity.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

#include "../Mat.hpp"

class GradientScheme
{
    public:

        GradientScheme() {}

        virtual ~GradientScheme() {}

        virtual void init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

        virtual Field<Mat<5,3>> calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const;
  
    protected:

        Field<Vars<3>> cellToCellDelta;

        void calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);
};

#endif // GRADIENTSCHEME_HPP