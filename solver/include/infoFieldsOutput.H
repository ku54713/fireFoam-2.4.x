    singleStepReactingMixture<gasHThermoPhysics>& singleMixture
    (
        dynamic_cast<singleStepReactingMixture<gasHThermoPhysics>&>
        (thermo)
    );

    // stoichiometric O2 to fuel ratio
    scalar s(singleMixture.s().value());

    // stoichiometric air to fuel ratio
    scalar stoicRatio(singleMixture.stoicRatio().value());  

    // heat of combustion [J/kg]
    scalar qFuel(singleMixture.qFuel().value());  

    label fuelIndex(singleMixture.fuelIndex());
    label inertIndex(singleMixture.inertIndex());

    const volScalarField* O2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
    if (thermo.composition().contains("O2")) 
    {
        O2Ptr = &thermo.composition().Y("O2");
    }
    const volScalarField& O2 = *O2Ptr;

    const volScalarField* CO2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
    if (thermo.composition().contains("CO2")) 
    {
        CO2Ptr = &thermo.composition().Y("CO2");
    }
    const volScalarField& CO2 = *CO2Ptr;

    const volScalarField& fu = Y[fuelIndex]; 

    scalar YO2Inf = 0.23301; //hardcode for now

    // Fuel molecular weight
    scalar Wu_ = singleMixture.speciesData()[fuelIndex].W();
    // Fuel enthalpy of formation
    scalar Hu_ = singleMixture.speciesData()[fuelIndex].hc();

    // compute stoichiometric mixture fraction
    scalar ftSt = 1.0 / ( 1.0 + stoicRatio );
    Info << "stoichiometric mixture fraction is = " << ftSt << nl;

    // create fileds for mixture fraction
    volScalarField ft
    (
     IOobject
     (
      "ft",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     (fu*s-O2+YO2Inf)/(s+YO2Inf)
    );

    // create fileds for patch integration of HRR
    surfaceScalarField HRR_fu
    ( 
        IOobject
        (
            "HRR_fu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(fu)*qFuel
    );

    // for outputing flame height
    OFstream outFlameHeight("outFlameHeight_"+runTime.timeName());
    outFlameHeight << "# Time    " << "FlameHeight  "  << nl;

    volScalarField flameHeight
    (
        IOobject
        (
            "flameHeight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );

    const pointField& cellCentres = mesh.cellCentres();


    // for check ft conservation
    surfaceScalarField phiFt
    ( 
        IOobject
        (
            "phiFt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiO2
    ( 
        IOobject
        (
            "phiO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiCO2
    ( 
        IOobject
        (
            "phiCO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    // create fileds for surface integration of sensible enthalpy flux
    surfaceScalarField phiHs
    (
        IOobject
        (
            "phiHs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(thermo.he()) 
    );

    // create fileds for surface integration of chemical enthalpy flux
    volScalarField hc(thermo.hc()); 
    surfaceScalarField phiHc
    (
        IOobject
        (
           "phiHc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create fileds for surface integration of total enthalpy flux
    volScalarField h(hc + thermo.he());
    surfaceScalarField phiH
    (
        IOobject
        (
            "phiH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create UT field for average purpose
    volVectorField UT
    (
        IOobject
        (
            "UT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U*T
    );

    // create rhoU field for average purpose
    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    // create rhoT field for average purpose
    volScalarField rhoT
    (
        IOobject
        (
            "rhoT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*T
    );

    // create rhoFu field for average purpose
    volScalarField rhoFu
    (
        IOobject
        (
            "rhoFu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*fu
    );

    // create B field for average purpose
    volSymmTensorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(U)
    );


    surfaceScalarField convectiveHeatFlux_L
    (
        IOobject
        (
            "convectiveHeatFlux_L",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    surfaceScalarField convectiveHeatFlux_T
    (
        IOobject
        (
            "convectiveHeatFlux_T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

/*
    volScalarField wallConvectiveHeatFlux
    (
        IOobject
        (
            "wallConvectiveHeatFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );
*/

    volScalarField QcWallFunction
    (
        IOobject
        (
            "QcWallFunction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    IOdictionary QcWallFunctionDict
    (
        IOobject
        (
            "QcWallFunctionDict",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            //IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar QcFlame =
    QcWallFunctionDict.lookupOrDefault<scalar>("QcFlame", 16000.0);
    
    scalar QcThreshold =
    QcWallFunctionDict.lookupOrDefault<scalar>("QcThreshold", 5000.0);

    volScalarField alphaM
    (
        IOobject
        (
            "alphaM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->alpha()
    );

    volScalarField muM
    (
        IOobject
        (
            "muM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->mu()
    );

    volScalarField d
    (
        IOobject
        (
            "d",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->alpha()/lewisNo/rho
    );

    volScalarField dSgs
    (
        IOobject
        (
            "dSgs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (turbulence->alphaEff()-turbulence->alpha())/rho
    );

