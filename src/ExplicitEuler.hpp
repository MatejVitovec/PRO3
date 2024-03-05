#ifndef EXPLICITEULER_HPP
#define EXPLICITEULER_HPP

#include "FVMScheme.hpp"

class ExplicitEuler : public FVMScheme
{
    public:

        //ExplicitEuler() : FVMScheme() {}
        ExplicitEuler(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : FVMScheme(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_)) {}

        virtual ~ExplicitEuler() {}

        void solve();        
        

    private:
    

};



#endif // EXPLICITEULER_HPP