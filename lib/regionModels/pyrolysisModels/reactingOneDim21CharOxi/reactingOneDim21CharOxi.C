/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "reactingOneDim21CharOxi.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "absorptionEmissionModel.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDim21CharOxi, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim21CharOxi, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim21CharOxi, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDim21CharOxi::updateCharOxi()
{

    const label charIndex = solidThermo_.composition().species()["char"];

    scalar rhoChar = solidThermo_.composition().rho(charIndex, 10000, 298);

    //Info << "mw of char is " <<  solidThermo_.composition().W(charIndex);

    // hardcode molecular weight for C and O2
    scalar mWO2 = 32;
    scalar mWChar = 12;
    scalar mWCO2 = 44;

    // hardcode HoC for char [J/kg]
    dimensionedScalar HocChar("HocChar",dimEnergy/dimMass,32.8e6);

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        scalarField& Xcharp = Xchar_.boundaryField()[patchI];
        scalarField& mCharp = mChar_.boundaryField()[patchI];
        scalarField& mCharBurntp = mCharBurnt_.boundaryField()[patchI];
        scalarField& charOxiShp = charOxiSh_.boundaryField()[patchI];
        scalarField& phiO2p = phiO2_.boundaryField()[patchI];
        scalarField& phiCO2p = phiCO2_.boundaryField()[patchI];

        const fvPatch& patch = regionMesh().boundary()[patchI];

        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            patch.patch()
        );
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];

        fvPatchScalarField O2 = nbrPatch.lookupPatchField<volScalarField, scalar>("O2");
        scalarField alpha = nbrPatch.lookupPatchField<volScalarField, scalar>("thermo:alpha");    //for now use molecular alpha
        scalarField O2Int(O2.patchInternalField());
        scalarField alphaDelta(alpha * nbrPatch.deltaCoeffs());

        mpp.distribute(O2);
        mpp.distribute(O2Int);
        mpp.distribute(alphaDelta);

        phiO2p = - alphaDelta * (O2Int - 0.0) * patch.magSf(); //[kg/s] negative
        phiCO2p = - phiO2p * mWCO2 / mWO2;  // positive, same as phiGas

        scalarField deltaMO2(- phiO2p * time_.deltaT().value());   //[kg]

        const scalarField& cellV = regionMesh().V();

        forAll(Xcharp, faceI)
        {
            const labelList& cells = boundaryFaceCells_[faceI];
            scalar charVolume = 0.0;
            scalar totalVolume = 0.0;
            forAll(cells, k)
            {
                const label cellI = cells[k];
                scalar X = Ys_[charIndex][cellI] * rho_[cellI] / rhoChar;

                charVolume += X * cellV[cellI];
                totalVolume += cellV[cellI];
            }
            Xcharp[faceI] = charVolume / totalVolume;
            mCharp[faceI] = rhoChar * charVolume;

            scalar charAvai = mCharp[faceI] - mCharBurntp[faceI];
            scalar dmCharBurnt = 0.0;
            if (charAvai < deltaMO2[faceI] / mWO2 * mWChar)
            {
                dmCharBurnt = charAvai;
                phiO2p[faceI] = - dmCharBurnt / mWChar * mWO2;
                phiCO2p[faceI] = dmCharBurnt / mWChar * mWCO2;
            }
            else
            {
                dmCharBurnt = deltaMO2[faceI] / mWO2 * mWChar;
            }
            mCharBurntp[faceI] += dmCharBurnt;

            //compute energy release rate in the first cell [kW/m3]
            charOxiSh_[cells[0]] = HocChar.value() * dmCharBurnt
                                   / cellV[cells[0]]
                                   / time_.deltaT().value();

            //for HRR diagnostics (sum at the patch is the HRR) [kW]
            charOxiShp[faceI] = charOxiSh_[cells[0]] * cellV[cells[0]];
        }
    }
}

void reactingOneDim21CharOxi::updateFields()
{
    reactingOneDim21::updateFields();

    updateCharOxi();
}

void reactingOneDim21CharOxi::solveEnergy()
{
    Info<< "reactingOneDim21CharOxi::solveEnergy()" << endl;

    if (debug)
    {
        Info<< "reactingOneDim21CharOxi::solveEnergy()" << endl;
    }

    tmp<volScalarField> alpha(solidThermo_.alpha());

    dimensionedScalar Cp0("Cp0", dimEnergy/dimMass/dimTemperature, solidThermo_.composition().Cp(0, 1.0e05, 300.) );
    dimensionedScalar Cp1("Cp1", dimEnergy/dimMass/dimTemperature, solidThermo_.composition().Cp(1, 1.0e05, 300.) );

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      - fvm::laplacian(alpha, h_)
      + fvc::laplacian(alpha, h_)
      - fvc::laplacian(kappa(), T())
     ==
        chemistrySh_
//      - fvm::Sp(solidChemistry_->RRg(), h_)
      + solidChemistry_->RRs(0)*T()*Cp0
      + solidChemistry_->RRs(1)*T()*Cp1
      + charOxiSh_
    );

    if (gasHSource_)
    {
        const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));
        hEqn += fvc::div(phiGas);
    }

    if (QrHSource_)
    {
        const surfaceScalarField phiQr(fvc::interpolate(Qr_)*nMagSf());
        hEqn += fvc::div(phiQr);
    }

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn += fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDim21CharOxi::reactingOneDim21CharOxi
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    reactingOneDim21(modelType, mesh, regionType),

    Xchar_
    (
        IOobject
        (
            "Xcharp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mChar_
    (
        IOobject
        (
            "mChar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    charOxiSh_
    (
        IOobject
        (
            "charOxiSh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    phiO2_
    (
        IOobject
        (
            "phiO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiCO2_
    (
        IOobject
        (
            "phiCO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    )

{}


reactingOneDim21CharOxi::reactingOneDim21CharOxi
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    reactingOneDim21(modelType, mesh, dict, regionType),

    Xchar_
    (
        IOobject
        (
            "Xcharp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mChar_
    (
        IOobject
        (
            "mChar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    charOxiSh_
    (
        IOobject
        (
            "charOxiSh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    phiO2_
    (
        IOobject
        (
            "phiO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiCO2_
    (
        IOobject
        (
            "phiCO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDim21CharOxi::~reactingOneDim21CharOxi()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
void reactingOneDim21CharOxi::info() const
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
