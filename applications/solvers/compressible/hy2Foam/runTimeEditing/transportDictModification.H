
if(thermo.transportDictionary().subDict("transportModels").lookupOrDefault<bool>("applyChanges", false))
{
    //---------------------------------------------------------
    // Re-initialisation of the transport quantities
    //---------------------------------------------------------
    transport.clear();
    transport = mixingRule::New
    (
        thermo,
        turbulence()
    );
    
    //---------------------------------------------------------
    // Re-definition of the multi-species transport model
    //---------------------------------------------------------
    speciesDiffusion.clear();
    speciesDiffusion = multiSpeciesTransportModel::New
    (
        thermo,
        turbulence()
    );
    
    speciesDiffusionName = thermo.transportDictionary().subDict("transportModels")
        .lookupOrDefault<word>("multiSpeciesTransport", word::null);
        
    writeWallHeatFlux = thermo.transportDictionary().subDict("transportModels").lookupOrDefault<bool>("writeWallHeatFlux", false);

    thermo.transportDictionary().subDict("transportModels").set<bool>("applyChanges", false);
    
    noSubRestart += 1;
}


if(thermo.transportDictionary().subDict("rarefiedParameters").lookupOrDefault<bool>("applyChanges", false))
{
    //---------------------------------------------------------
    // Re-definition of the rarefaction parameters
    //---------------------------------------------------------
    rarefactionParameters.clear();
    rarefactionParameters = rarefactionParameter::New
    (
        thermo,
        turbulence()
    );
    
    rarefactionParameters->correct(U);
    
    thermo.transportDictionary().subDict("rarefiedParameters").set<bool>("applyChanges", false);
    
    noSubRestart += 1;
}

