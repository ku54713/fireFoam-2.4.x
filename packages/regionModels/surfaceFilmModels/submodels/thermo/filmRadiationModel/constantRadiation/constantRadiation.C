/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "constantRadiation.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    constantRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantRadiation::constantRadiation
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, owner, dict),
    QrConst_
    (
        IOobject
        (
            typeName + "_QrConst",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh()
    ),
    mask_
    (
        IOobject
        (
            typeName + "_mask",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("one", dimless, 1.0)//,
        // zeroGradientFvPatchScalarField::typeName
    ),
    absorptivity_(readScalar(coeffDict_.lookup("absorptivity"))),
    Qin_
    (
        IOobject
        (
            "Qin", // same name as Qin on primary region to enable mapping
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
        zeroGradientFvPatchScalarField::typeName
        ),
    timeStart_(readScalar(coeffDict_.lookup("timeStart"))),
    duration_(readScalar(coeffDict_.lookup("duration")))
{
    mask_ = pos(mask_ - 0.5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantRadiation::~constantRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void constantRadiation::correct()
{}


tmp<volScalarField> constantRadiation::Shs()
{
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "_Shs",
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    const scalar time = owner().time().value();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        scalarField& Shs = tShs();
        const scalarField& Qr = QrConst_.internalField();
        // const scalarField& alpha = owner_.alpha().internalField();

        // Shs = mask_*Qr*alpha*absorptivity_;
        Shs = mask_*Qr*absorptivity_;
        // Qin_ is used in turbulentTemperatureRadiationCoupledMixedST 
        // boundary condition
        Qin_.internalField() = mask_*Qr;
        Qin_.correctBoundaryConditions();
    }

    return tShs;
}

tmp<volScalarField> constantRadiation::ShsConst() const
{
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "_Shs",
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    const scalar time = owner().time().value();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        scalarField& Shs = tShs();
        const scalarField& Qr = QrConst_.internalField();
        const scalarField& alpha = owner_.alpha().internalField();

        // Shs = mask_*Qr*alpha*absorptivity_;
        Shs = mask_*Qr*absorptivity_;
        // Qin_ is used in turbulentTemperatureRadiationCoupledMixedST 
        // boundary condition
    }

    return tShs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
