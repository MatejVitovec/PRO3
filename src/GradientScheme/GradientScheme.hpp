#ifndef GRADIENTSCHEME_HPP
#define GRADIENTSCHEME_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"

class GradientScheme
{
    public:

        GradientScheme() {}

        virtual void init(const Mesh& mesh);

        virtual Field<std::array<Vars<5>, 3>> calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const;


    protected:

};

#endif // GRADIENTSCHEME_HPP