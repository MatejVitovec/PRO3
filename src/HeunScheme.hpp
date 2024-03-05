#ifndef HEUNSCHEME_HPP
#define HEUNSCHEME_HPP

#include "FVMScheme.hpp"

class HeunScheme : public FVMScheme
{
    public:

        //HeunScheme() : FVMScheme() {}
        HeunScheme(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : FVMScheme(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_)) {}


        virtual ~HeunScheme() {}

        void solve();        
        

    private:
     

};



#endif // HEUNSCHEME_HPP