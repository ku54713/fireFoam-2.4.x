    volScalarField PsiS
    (
        IOobject
        (
            "PsiS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("PsiS", dimless, scalar(0))
    );
    volScalarField qrflamelet
    (
        IOobject
        (
            "qrflamelet",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("qrflamelet", dimMass/(dimLength*pow3(dimTime)), scalar(0))
    );
    volScalarField chiSoot
    (
        IOobject
        (
            "chiSoot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("chiSoot", dimless, scalar(0))
    );
    volScalarField PsiR
    (
        IOobject
        (
            "PsiR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("PsiR", dimless, scalar(0))
    );
    volScalarField strainR
    (
        IOobject
        (
            "strainR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("strainR", pow(dimTime,-1), scalar(0))
    );
    volScalarField srC
    (
        IOobject
        (
            "srC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("srC", pow(dimTime,-1), scalar(0))
    );
    dimensionedScalar C_src
    (
      	"C_src",
      	pow(dimTime,-1),
      	1.0
    );
    dimensionedScalar C_eta
    (
      	"C_eta",
      	dimless,
      	1.0
    );
    volScalarField Hp
    (
        IOobject
        (
            "Hp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("Hp", dimless, scalar(0))
    );


// soot dictionary
IOdictionary sootModelPropertiesDict
(
    IOobject
    (
        "sootModelProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
// soot dictionary default values
Switch solveSoot
(
    sootModelPropertiesDict.lookupOrDefault<Switch>("solveSoot",false)
);
if( solveSoot ) {
	Info << "Initializing soot model" << nl;
}
Switch constSR
(
    sootModelPropertiesDict.lookupOrDefault<Switch>("constSR",false)
);
Switch constHP
(
    sootModelPropertiesDict.lookupOrDefault<Switch>("constHP",false)
);
Switch newHP
(
    sootModelPropertiesDict.lookupOrDefault<Switch>("newHP",true)
);

Switch radScaling
(
    sootModelPropertiesDict.lookupOrDefault<Switch>("radScaling",false)
);

wordList radSc_fPatchName;
wordList radSc_pPatchName;
if( radScaling ) {
const wordList radSc_fPatchName_
(
	sootModelPropertiesDict.lookup("radSc_fPatch")
);
const wordList radSc_pPatchName_
(
	sootModelPropertiesDict.lookup("radSc_pPatch")
);
radSc_fPatchName = radSc_fPatchName_;
radSc_pPatchName = radSc_pPatchName_;
}

scalar radSc_Xr_f
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("radSc_Xr_f",0.3)
);
scalar radSc_Xr_p
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("radSc_Xr_p",0.2)
);
scalar radSc_Ys_f
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("radSc_Ys_f",0.02)
);
scalar radSc_Ys_p
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("radSc_Ys_p",0.01)
);
scalar Csc
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("Csc",1)
);
scalar Ceta
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("Ceta",0.04)
);
scalar cSR
(
    sootModelPropertiesDict.lookupOrDefault<scalar>("cSR",10.0)
);

// log messages
if( solveSoot ) {
if( constSR ) {
	Info << "	-> constant strain rate = " << cSR << " 1/s" << endl;
}
else {
	Info << "	-> Ceta = " << Ceta << endl;
}
if( constHP ) {
	Info << "	-> zero prior heat loss" << nl;
}
if( radScaling ) {
	Info << "	-> multi-fuel radiation scaling" << nl;
	Info << "	   (Xr,fuel = " << radSc_Xr_f << ", Ys,fuel = " << radSc_Ys_f << ", Xr,pyro = " << radSc_Xr_p << ", Ys,pyro = " << radSc_Ys_p << ")" << nl;
}
Info << "" << nl;
}

// lookup table
#include "sootLookup.H"
