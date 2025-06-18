
typedef void (*FnPtr)
(
    const bool&,
    const bool&,
    // const bool&,
    // const scalar&,
    const vector&,
    const scalar&,  // new
    const scalar&,  // new
    const scalar&,  // new
    const scalar&,  // new
    const scalar&,  // new
    const scalar&,
    const scalar&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const vector&,
    const vector&,
    const vector&,
    const vector&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const scalar&,
    const scalar&,
    scalar&,
    vector&,
    scalar&,
    vector&,
    PtrList<scalar>&,
    scalar&
);

std::map<word, FnPtr> fluxMap;
fluxMap["Kurganov"] = KurganovFlux;
fluxMap["Tadmor"] = TadmorFlux;
fluxMap["KFVS"] = KFVSFlux;
fluxMap["LxF"] = LxFFlux;
fluxMap["ROE"] = ROEFlux;
fluxMap["ROECorrect"] = ROECorrectFlux;
fluxMap["HLLC"] = HLLCFlux;
fluxMap["AUSMPLUS"] = AUSMPLUSFlux;
fluxMap["AUSMPLUSUP"] = AUSMPLUSUPFlux;
fluxMap["AUSMPLUSUP0"] = AUSMPLUSUP0Flux;

word fluxScheme("Kurganov");
FixedList<word, 2> blendedSchemes;

if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if
    (
         (fluxScheme == "Tadmor") 
        or (fluxScheme == "Kurganov")
        or (fluxScheme == "AUSMPLUS")
        or (fluxScheme == "AUSMPLUSUP")
        or (fluxScheme == "AUSMPLUSUP0")
        or (fluxScheme == "KFVS")
        or (fluxScheme == "LxF")
        or (fluxScheme == "ROE")
        or (fluxScheme == "ROECorrect")
        or (fluxScheme == "HLLC")
    )
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else if (fluxScheme == "blended")
    {
        if (mesh.schemesDict().readIfPresent("blendedSchemes", blendedSchemes))
        {
            Info<< "fluxScheme for shock vicinity: " << blendedSchemes[0]
                << endl;
            Info<< "fluxScheme for near-wall region: " << blendedSchemes[1]
                << endl;
        }
        else
        {
            FatalErrorIn
            (
                "readFluxScheme"
            )   << "blendedSchemes: " << fluxScheme
                << " is not a valid choice. "
                << "Options are: Kurganov, Tadmor, ROECorrect, AUSMPLUS, AUSMPLUSUP and AUSMPLUSUP0 (ROE, KFVS, LxF, HLLC)."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "readFluxScheme"
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. Options are: Kurganov, Tadmor, ROECorrect, AUSMPLUS, AUSMPLUSUP, AUSMPLUSUP0 (ROE, KFVS, LxF, HLLC), "
            << "and blended"
            << abort(FatalError);
    }
}

// scalar sqrMachInf =
//     sqr(mesh.schemesDict().lookupOrDefault<scalar>("Mach_inf", SMALL));


// --- Declarations ---
surfaceScalarField amaxSf
(
    IOobject
    (
        "amaxSf",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("amaxSf", dimVolume/dimTime, 0.0)
);

surfaceVectorField aByU
(
    IOobject
    (
        "aByU",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("aByU", dimVelocity, vector::zero)
);

surfaceVectorField phiUp
(
    IOobject
    (
        "phiUp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("phiUp", dimMass/dimTime*dimVelocity, vector::zero)
);

surfaceScalarField phiEp
(
    IOobject
    (
        "phiEp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("phiEp", dimEnergy/dimTime, 0.0)
);

PtrList<surfaceScalarField> phiEvk(noSolvedVibEq);

if (downgradeSingleT)
{
    phiEvk.clear();
}
else if (downgradeSingleTv)
{
    phiEvk.resize(1);

    phiEvk.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                "phiEvk",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("phiEvk", dimEnergy/dimTime, 0.0)
        )
    );
}
else
{
    forAll(phiEvk, speciei)
    {
        phiEvk.set
        (
            speciei,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiEvk" + Y[speciei].name(),
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "phiEvk" + Y[speciei].name(),
                    dimEnergy/dimTime,
                    0.0
                )
            )
        );
    }
}

