#include "Hllc.hpp"

Vars<5> Hllc::claculateFlux(const Compressible& wl, const Compressible& wr, const Vars<3>& normalVector) const
{    
    double rhoL = wl.density();
    double pL = wl.pressure();
    double nuL = wl.normalVelocity(normalVector);
    double aL = wl.soundSpeed();

    double rhoR = wr.density();
    double pR = wr.pressure();
    double nuR = wr.normalVelocity(normalVector);
    double aR = wr.soundSpeed();
    
    //PVRS
    double pm = std::fmax(0, 0.5*(pL + pR) - 0.5*(nuR - nuL)*0.5*(rhoL + rhoR)*0.5*(aL + aR));

    double sl;
    double sr;
    double sm;

    if (pm <= pL)
    {
        sl = nuL - aL;
    }        
    else
    {
        sl = 0.5*(nuL + nuR) - 0.5*(aL + aR);
    }
        
    if (pm <= pR)
    {
        sr = nuR + aR;
    }        
    else
    {
        sr = 0.5*(nuL + nuR) + 0.5*(aL + aR);
    }
        
    //contact wave speed
    sm = pR - pL + rhoL*nuL*(sl - nuL) - rhoR*nuR*(sr - nuR);
    sm = sm / (rhoL*(sl - nuL) - rhoR*(sr - nuR));


    //HLLC scheme
    if (sl >= 0)
    {
        //left state
        return wl.flux(normalVector);
    }
    else if (sr <= 0)
    {
        //right state
        return wr.flux(normalVector);
    }
    else if (sm >= 0)
    {
        //middle-left state

        double uL = wl.velocityU();
        double vL = wl.velocityV();
        double wL = wl.velocityW();
        double EL = wl.totalEnergy();

        pm = pL + rhoL*(sl - nuL)*(sm - nuL);
        double rhoM = rhoL*(sl - nuL)/(sl - sm);

        return Vars<5>({rhoL*nuL + sl*(rhoM - rhoL),
                        rhoL*nuL*uL + pL*normalVector[0] + sl*((rhoM - rhoL)*uL + (pm - pL)/(sl - sm)*normalVector[0]),
                        rhoL*nuL*vL + pL*normalVector[1] + sl*((rhoM - rhoL)*vL + (pm - pL)/(sl - sm)*normalVector[1]),
                        rhoL*nuL*wL + pL*normalVector[2] + sl*((rhoM - rhoL)*wL + (pm - pL)/(sl - sm)*normalVector[2]),
                        rhoL*nuL*(EL + pL/rhoL) + sl*((rhoM - rhoL)*EL + (pm*sm - pL*nuL)/(sl - sm))});
    }
    else
    {
        //middle-right state

        double uR = wr.velocityU();
        double vR = wr.velocityV();
        double wR = wr.velocityW();
        double ER = wr.totalEnergy();

        pm = pR + rhoR*(sr - nuR)*(sm - nuR);
        double rhoM = rhoR*(sr - nuR)/(sr - sm);
        return Vars<5>({rhoR*nuR + sr*(rhoM - rhoR),
                        rhoR*nuR*uR + pR*normalVector[0] + sr*((rhoM - rhoR)*uR + (pm - pR)/(sr - sm)*normalVector[0]),
                        rhoR*nuR*vR + pR*normalVector[1] + sr*((rhoM - rhoR)*vR + (pm - pR)/(sr - sm)*normalVector[1]),
                        rhoR*nuR*wR + pR*normalVector[2] + sr*((rhoM - rhoR)*wR + (pm - pR)/(sr - sm)*normalVector[2]),
                        rhoR*nuR*(ER + pR/rhoR) + sr*((rhoM - rhoR)*ER + (pm*sm - pR*nuR)/(sr - sm))});
    }
}