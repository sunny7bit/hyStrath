// --- Directed interpolation of primitive fields onto faces
tmp<surfaceScalarField> rho_pos_tmp, rho_neg_tmp;
fvc::surfaceReconstruct(rho, rho_pos_tmp, rho_neg_tmp, "reconstruct(rho)");
surfaceScalarField rho_pos = rho_pos_tmp();  
surfaceScalarField rho_neg = rho_neg_tmp();  

PtrList<surfaceScalarField> pD_pos(noSolvedVibEq);
PtrList<surfaceScalarField> pD_neg(noSolvedVibEq);

if(downgradeSingleT or downgradeSingleTv)
{
    pD_pos.clear();
    pD_neg.clear();
}
else
{
    forAll(pD_pos, speciei)
    { 
        pD_pos.set
        (
            speciei,
            new surfaceScalarField
            (
                "pD_pos_" + Y[speciei].name(),
                fvc::surfaceReconstruct(pD[speciei], pos, "reconstruct(rho)")
            )
        );
    }

    forAll(pD_neg, speciei)
    { 
        pD_neg.set
        (
            speciei,
            new surfaceScalarField
            (
                "pD_neg_" + Y[speciei].name(),
                fvc::surfaceReconstruct(pD[speciei], neg, "reconstruct(rho)")
            )
        );
    }
}

// ---
tmp<surfaceVectorField> U_pos_tmp, U_neg_tmp;
fvc::surfaceReconstruct(U, U_pos_tmp, U_neg_tmp, "reconstruct(U)");
surfaceVectorField U_pos = U_pos_tmp();  
surfaceVectorField U_neg = U_neg_tmp();  
// surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
// surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

tmp<surfaceScalarField> p_pos_tmp, p_neg_tmp;
fvc::surfaceReconstruct(p, p_pos_tmp, p_neg_tmp, "reconstruct(rho)");
surfaceScalarField p_pos = p_pos_tmp();  
surfaceScalarField p_neg = p_neg_tmp();  

// surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
// surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

// tmp<surfaceVectorField> rhoU_pos_tmp, rhoU_neg_tmp;
// fvc::surfaceReconstruct(rhoU, rhoU_pos_tmp, rhoU_neg_tmp, "reconstruct(U)");
//  surfaceVectorField& rhoU_pos = rhoU_pos_tmp(); 
//  surfaceVectorField& rhoU_neg = rhoU_neg_tmp();
surfaceVectorField rhoU_pos("rhoU_pos", rho_pos*U_pos);
surfaceVectorField rhoU_neg("rhoU_neg", rho_neg*U_neg);

volScalarField rPsi(1.0/psi);
// tmp<surfaceScalarField> rPsi_pos_tmp, rPsi_neg_tmp;
// fvc::surfaceReconstruct(rPsi, rPsi_pos_tmp, rPsi_neg_tmp, "reconstruct(T)");
//  surfaceScalarField& rPsi_pos = rPsi_pos_tmp(); 
//  surfaceScalarField& rPsi_neg = rPsi_neg_tmp(); 
surfaceScalarField rPsi_pos("rPsi_pos", p_pos/rho_pos);
surfaceScalarField rPsi_neg("rPsi_neg", p_neg/rho_neg);

tmp<surfaceScalarField> e_pos_tmp, e_neg_tmp;
fvc::surfaceReconstruct(e, e_pos_tmp, e_neg_tmp, "reconstruct(T)");
surfaceScalarField e_pos = e_pos_tmp(); 
surfaceScalarField e_neg = e_neg_tmp(); 

PtrList<surfaceScalarField> evk_pos(noSolvedVibEq);
PtrList<surfaceScalarField> evk_neg(noSolvedVibEq);

if(downgradeSingleT)
{
    evk_pos.clear();
    evk_neg.clear();  
}
else if(downgradeSingleTv)
{
    evk_pos.resize(1);
    evk_neg.resize(1);
    
    evk_pos.set
    (
        0,
        new surfaceScalarField
        (
            "ev_pos",
            fvc::surfaceReconstruct(ev, pos, "reconstruct(T)")
        )
    );
    
    evk_neg.set
    (
        0,
        new surfaceScalarField
        (
            "ev_neg",
            fvc::surfaceReconstruct(ev, neg, "reconstruct(T)")
        )
    );
}
else
{
    forAll(evk_pos, speciei)
    { 
        evk_pos.set
        (
            speciei,
            new surfaceScalarField
            (
                "evk_pos_" + Y[speciei].name(),
                fvc::surfaceReconstruct(hevk[speciei], pos, "reconstruct(T)")
            )
        );
    }

    forAll(evk_neg, speciei)
    { 
        evk_neg.set
        (
            speciei,
            new surfaceScalarField
            (
                "evk_neg_" + Y[speciei].name(),
                fvc::surfaceReconstruct(hevk[speciei], neg, "reconstruct(T)")
            )
        );
    }
}

volScalarField gammatr("gammatr", thermo.Cp_t()/thermo.Cv_t()); // NEW VINCENT (also useful for bdry cdts)
volScalarField c(sqrt(gammatr*rPsi)); // NEW VINCENT (see Candler 2009)
// volScalarField c(sqrt(gammatr*p/rho)); // NEW VINCENT (see Candler 2009)
Mach = mag(U)/c;

tmp<surfaceScalarField> c_pos_tmp, c_neg_tmp;
fvc::surfaceReconstruct(c, c_pos_tmp, c_neg_tmp, "reconstruct(T)");
surfaceScalarField c_pos = c_pos_tmp(); 
surfaceScalarField c_neg = c_neg_tmp(); 

surfaceScalarField cSf_pos //- Intermediate quantity to appear in Eq.8a
(
    "cSf_pos",
    c_pos*mesh.magSf()
);
// cSf_pos.setOriented();
surfaceScalarField cSf_neg //- Intermediate quantity to appear in Eq.8b
(
    "cSf_neg",
    c_neg*mesh.magSf()
);
// cSf_neg.setOriented();

tmp<surfaceScalarField> gamma_pos_tmp, gamma_neg_tmp;
fvc::surfaceReconstruct(gammatr, gamma_pos_tmp, gamma_neg_tmp, "reconstruct(T)");
surfaceScalarField gamma_pos = gamma_pos_tmp(); 
surfaceScalarField gamma_neg = gamma_neg_tmp(); 