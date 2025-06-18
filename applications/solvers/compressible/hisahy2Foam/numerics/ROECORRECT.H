/*---------------------------------------------------------------------------*\
Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix. This is a first order scheme
\*---------------------------------------------------------------------------*/
void ROECorrectFlux 
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
    
    // const scalar sum_evk_pos = 
    //     std::accumulate(evk_pos.begin(), evk_pos.end(), 0.0);
        
    // const scalar etr_pos = e_pos - sum_evk_pos;

    // const scalar gammaMinusOne = p_pos/(rho_pos*etr_pos);

    //- Enthalpy
    const scalar h_pos  = e_pos + p_pos/max(rho_pos,SMALL) + 0.5 * magSqr(U_pos);
    const scalar h_neg  = e_neg + p_neg/max(rho_neg,SMALL) + 0.5 * magSqr(U_neg);

    const scalar rho_pos_sqrt = Foam::sqrt(max(rho_pos,SMALL));
    const scalar rho_neg_sqrt = Foam::sqrt(max(rho_neg,SMALL));
    const scalar fact_pos = rho_pos_sqrt/(rho_pos_sqrt + rho_neg_sqrt);
    const scalar fact_neg = 1.0 - fact_pos;

    //- Roe average state
    const scalar rhoTilde  = rho_pos_sqrt*rho_neg_sqrt;
    const vector UTilde = U_pos*fact_pos + U_neg*fact_neg;
    const scalar hTilde = h_pos*fact_pos + h_neg*fact_neg;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar gammaTilde = gamma_pos*fact_pos + gamma_neg*fact_neg;
 
    const scalar contrVTilde = UTilde & normalVector;
    const scalar cTilde = Foam::sqrt(max((gammaTilde -1)*(hTilde - 0.5*qTildeSquare),SMALL)); // TODO
    // neg---right;left---pos
    const scalar deltaP = p_neg - p_pos;
    const scalar vel_pos_normal = U_pos & normalVector;
    const scalar vel_neg_normal = U_neg & normalVector;
    const vector deltaU = U_neg - U_pos;
    const scalar deltaContrV = vel_neg_normal - vel_pos_normal;
    const scalar deltaRho = rho_neg - rho_pos;
    
    // Roe and Pike - formulation
    const scalar r1 =
        (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
    const scalar r2 = deltaRho - deltaP/sqr(cTilde);
    const scalar r3 =
        (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // rho row:
    const scalar l1rho = 1;
    const scalar l2rho = 1;
    const scalar l3rho = 0;
    const scalar l4rho = 1;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;
    // second U column
    const vector l2U = UTilde;
    // third U column
    const vector l3U = deltaU - deltaContrV*normalVector;
    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // derived from algebra by hand, only for Euler equation usefull
    scalar lambda1 = mag(contrVTilde - cTilde);
    scalar lambda2 = mag(contrVTilde);
    scalar lambda3 = mag(contrVTilde + cTilde);

    // entropy correction: Felipe Portela, 9/Oct/2013
    const scalar cpos = Foam::sqrt
    (
        max
        (
            (gamma_pos - 1)*(h_pos - 0.5*magSqr(U_pos)),
            SMALL
        )
    );
    const scalar cneg = Foam::sqrt
    (
        max
        (
            (gamma_neg - 1)*(h_neg - 0.5*magSqr(U_neg)),
            SMALL
        )
    );

    // First eigenvalue: U - c
    scalar eps = 2*max(0,(vel_neg_normal - cneg) - (vel_pos_normal - cpos));
    if (lambda1 < eps)
    {
        lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    }

    // Second eigenvalue: U
    eps = 2*max(0, vel_neg_normal - vel_pos_normal);
    if (lambda2 < eps)
    {
        lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    }

    // Third eigenvalue: U + c
    eps = 2*max(0,(vel_neg_normal + cneg) - (vel_pos_normal + cpos));
    if (lambda3 < eps)
    {
        lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    }

    // Components of deltaF1
    const scalar diffF11 = lambda1*r1*l1rho;
    const vector diffF124 = lambda1*r1*l1U;
    const scalar diffF15 = lambda1*r1*l1e;

    // Components of deltaF2
    const scalar diffF21 = lambda2*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambda2*(r2*l2e + rhoTilde*l3e);

    // Components of deltaF3
    const scalar diffF31 = lambda3*r3*l4rho;
    const vector diffF324 = lambda3*r3*l4U;
    const scalar diffF35 = lambda3*r3*l4e;

    // Left flux 5-vector
    const scalar fluxpos11 = rho_pos*vel_pos_normal;
    const vector fluxpos124 = U_pos*fluxpos11 + normalVector*p_pos;
    const scalar fluxpos15 = h_pos*fluxpos11;
    // Right flux 5-vector
    const scalar fluxneg11 = rho_neg*vel_neg_normal;
    const vector fluxneg124 = U_neg*fluxneg11 + normalVector*p_neg;
    const scalar fluxneg15 = h_neg*fluxneg11;

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxpos11 + fluxneg11 - (diffF11 + diffF21 + diffF31));

    const vector flux24 =
        0.5*(fluxpos124 + fluxneg124 - (diffF124 + diffF224 + diffF324));

    const scalar flux5 =
        0.5*(fluxpos15 + fluxneg15 - (diffF15 + diffF25 + diffF35));

    // Compute private data
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;

    // ---
    if(downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if(downgradeSingleTv)
    {
        rhoEvkFlux[0] = 0.5*(evk_pos[0]*fluxpos11 + evk_neg[0]*fluxneg11 - (diffF15 + diffF25 + diffF35))*magSf;
    }
    else
    {
        forAll(rhoEvkFlux, speciei)
        { 
            rhoEvkFlux[speciei] = 0.5*(evk_pos[speciei]*pD_pos[speciei]*vel_pos_normal + evk_neg[speciei]*pD_neg[speciei]*vel_neg_normal - (diffF15 + diffF25 + diffF35))*magSf;
        }
    }
}
