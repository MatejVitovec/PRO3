#ifndef EXPLICITEULER_HPP
#define EXPLICITEULER_HPP

#include "FVMScheme.hpp"

class ExplicitEuler : public FVMScheme
{
    public:

        //ExplicitEuler() : FVMScheme() {}
        ExplicitEuler(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_) : FVMScheme(std::move(mesh_), std::move(fluxSolver_)) {}


        virtual ~ExplicitEuler() {}

        void solve();        
        

    private:
    
        Field<Vars<5>> calculateResidual();
        Field<Compressible> explicitIntegration(const Field<Vars<5>>& res);
        

};



#endif // EXPLICITEULER_HPP