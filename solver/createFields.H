    Info<< "Creating combustion model\n" << nl;

    autoPtr<combustionModels::psiCombustionModel> combustion
    (
        combustionModels::psiCombustionModel::New
        (
            mesh
        )
    );

    Info<< "Reading thermophysical properties\n" << nl;

    //hsCombustionThermo& thermo = combustion->thermo();
    psiReactionThermo& thermo = combustion->thermo();

    SLGThermo slgThermo(mesh, thermo);

    basicMultiComponentMixture& composition = thermo.composition();

    PtrList<volScalarField>& Y = composition.Y();

    const word inertSpecie(thermo.lookup("inertSpecie"));

    #include "readAdditionalThermo.H"

    Info<< "Creating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    //volScalarField& hs = thermo.hs();
    volScalarField& hs = thermo.he();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    Info<< "Reading field p_rgh\n" << endl;
    
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading field U\n" << endl;
    
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, phi, thermo)
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

/*
IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
*/

volScalarField dQ
(
    IOobject
    (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimMass/pow3(dimTime)/dimLength, 0.0)
);


#include "infoFieldsOutput.H"


Info<< "Creating field DpDt\n" << nl;
volScalarField DpDt =
   fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


Info<< "Calculating field g.h\n" << nl;
volScalarField gh("gh", g & mesh.C());

    surfaceScalarField ghf("ghf", g & mesh.Cf());

    p = p_rgh + rho*gh;
    thermo.correct();
    rho = thermo.rho();

    // Force p_rgh to be consistent with p
    p_rgh = p - rho*gh;

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
    
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(hs);

    IOdictionary additionalControlsDict
    (
        IOobject
        (
            "additionalControls",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Switch solvePrimaryRegion
    (
        additionalControlsDict.lookup("solvePrimaryRegion")
    );
    Switch solvePyrolysisRegion
    (
        additionalControlsDict.lookupOrDefault<bool>("solvePyrolysisRegion", true)
    );


    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.Cp()*thermo.alpha()
    );
    
    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.Cp()
    );
    
