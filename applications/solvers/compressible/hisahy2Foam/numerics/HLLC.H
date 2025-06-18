#include <numeric>

/*---------------------------------------------------------------------------*\
Application
    hllcFlux

Description
    HLLC flux class, computes numerical fluxes based on HLLC scheme. Uses
    primitive variables as input and gives back conservative numerical fluxes.

    Toro, E.I. and Spruce, M. and Spears, W.: Restoration of the Contact
    Surface in the HLL-Riemann Solver, Shock Waves, 4:25-34, 1996
\*---------------------------------------------------------------------------*/

void HLLCFlux 
(
    const bool& downgradeSingleT,
    const bool& downgradeSingleTv,
    const vector& Sf,
    const scalar& magSf,
    const scalar& c_pos,
    const scalar& c_neg,
    const scalar& gamma_pos,
    const scalar& gamma_neg,
    const scalar& cSf_pos,
    const scalar& cSf_neg,
    const scalar& rho_pos,
    const scalar& rho_neg,
    const PtrList<scalar>& pD_pos,
    const PtrList<scalar>& pD_neg,
    const vector& U_pos,
    const vector& U_neg,
    const vector& rhoU_pos,
    const vector& rhoU_neg,
    const scalar& p_pos,
    const scalar& p_neg,
    const PtrList<scalar>& evk_pos,
    const PtrList<scalar>& evk_neg,
    const scalar& e_pos,
    const scalar& e_neg,
    scalar& amaxSf,
    vector& aByU,
    scalar& rhoFlux,
    vector& rhoUFlux,
    PtrList<scalar>& rhoEvkFlux,
    scalar& rhoEFlux
)
{
    // Maximum eigenvalue
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf)
        + 0.5*(cSf_pos + cSf_neg);
        
    aByU = (U_pos + U_neg)/2.0;    
        
    const vector normalVector = Sf/magSf;
    const scalar qpos  = U_pos & normalVector;
    const scalar qneg = U_neg & normalVector;

    // DensityTotalEnergy
    const scalar rhoE_pos = rho_pos*(e_pos + 0.5 * magSqr(U_pos));
    const scalar rhoE_neg = rho_neg*(e_neg + 0.5 * magSqr(U_neg));

    //- Enthalpy
    const scalar h_pos  = (rhoE_pos + p_pos)/rho_pos;
    const scalar h_neg  = (rhoE_neg + p_neg)/rho_neg;

    const scalar rho_pos_sqrt = Foam::sqrt(max(rho_pos,0.0));
    const scalar rho_neg_sqrt = Foam::sqrt(max(rho_neg,0.0));
    const scalar fact_pos = rho_pos_sqrt/stabilise((rho_pos_sqrt + rho_neg_sqrt),VSMALL);
    const scalar fact_neg = 1.0 - fact_pos;

    //- Roe average state
    const vector UTilde = U_pos*fact_pos + U_neg*fact_neg;
    const scalar hTilde = h_pos*fact_pos + h_neg*fact_neg;
    const scalar gammaTilde = gamma_pos*fact_pos + gamma_neg*fact_neg;
 
    const scalar contrUTilde = UTilde & normalVector;
    const scalar cTilde = Foam::sqrt(max((gammaTilde -1)*(hTilde - 0.5*sqr(contrUTilde)),0.0)); // TODO

    // Step 3: compute signal speeds for face:
    const scalar Spos  = min(qpos-c_pos, contrUTilde-cTilde);
    const scalar Sneg = max(contrUTilde+cTilde, qneg+c_neg);

    const scalar SStar = (rho_neg*qneg*(Sneg-qneg)
      - rho_pos*qpos*(Spos - qpos) + p_pos - p_neg)/
        stabilise((rho_neg*(Sneg-qneg)-rho_pos*(Spos - qpos)),VSMALL);

    // Compute pressure in star region from the right side
    const scalar pStarneg =
        rho_neg*(qneg - Sneg)*(qneg - SStar) + p_neg;

    // Should be equal to the left side
    const scalar pStarpos  =
        rho_pos*(qpos -  Spos)*(qpos - SStar) + p_pos;
    // Give a warning if this is not the case
    if (mag(pStarneg - pStarpos) > 1e-6)
    {
        Info << "mag(pStarneg - pStarpos) > VSMALL " << endl;
    }
    // Use pStarneg for pStar, as in theory, pStarneg == pStarpos
    const scalar pStar = pStarneg;
    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;
    if (Spos>=0)
    {
        // compute F_l
        convectionSpeed = qpos;
        rhoState  = rho_pos;
        rhoUState = rhoU_pos;
        rhoEState = rhoE_pos;
        pState = p_pos;
    }
    else if (SStar>=0)
    {
        scalar omegapos = scalar(1.0)/stabilise((Spos - SStar), VSMALL);

        // Compute left star region
        convectionSpeed = SStar;
        rhoState  = omegapos*(Spos - qpos)*rho_pos;
        rhoUState = omegapos*((Spos - qpos)*rhoU_pos
        + (pStar - p_pos)*normalVector);
        rhoEState = omegapos*((Spos - qpos)*rhoE_pos
        - p_pos*qpos + pStar*SStar);
        pState = pStar;
    }
    else if (Sneg>=0)
    {
        scalar omeganeg = scalar(1.0)/stabilise((Sneg - SStar), VSMALL);

        // compute right star region
        convectionSpeed = SStar;
        rhoState  = omeganeg*(Sneg - qneg)*rho_neg;
        rhoUState = omeganeg*((Sneg - qneg)*rhoU_neg
        + (pStar - p_neg)*normalVector);
        rhoEState = omeganeg*((Sneg - qneg)*rhoE_neg
        - p_neg*qneg + pStar*SStar);
        pState = pStar;
    }
    else if (Sneg<0)
    {
        // compute F_r
        convectionSpeed = qneg;
        rhoState  = rho_neg;
        rhoUState = rhoU_neg;
        rhoEState = rhoE_neg;
        pState = p_neg;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState))*magSf;
}
