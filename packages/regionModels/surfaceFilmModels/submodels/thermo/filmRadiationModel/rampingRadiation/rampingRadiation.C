/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "rampingRadiation.H"
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

defineTypeNameAndDebug(rampingRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    rampingRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rampingRadiation::rampingRadiation
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh()//,
        // dimensionedScalar("one", dimless, 1.0)//,
        // zeroGradientFvPatchScalarField::typeName
    ),
    absorptivity_(readScalar(coeffDict_.lookup("absorptivity"))),
    // initialValue_(readScalar(coeffDict_.lookup("initialValue"))),
    rampTimeInterval_(readScalar(coeffDict_.lookup("rampTimeInterval"))),
    rampStep_(readScalar(coeffDict_.lookup("rampStep"))),
    rampStartTime_(0.0),
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
    rampStartTime_=max(timeStart_,owner.time().value());

    // for non-uniform heat flux (2d parabolic profile)
    // forAll(QrConst_,cellI){
    //     scalar x=owner.regionMesh().cellCentres()[cellI][0];
    //     scalar y=owner.regionMesh().cellCentres()[cellI][1];
    //     scalar xCenter = 0.255;
    //     scalar yCenter = 0.605;
    //     QrConst_[cellI] = (-26.02*pow(x-xCenter,2)-25.94*pow(y-yCenter,2)+54.98)*1e3;
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampingRadiation::~rampingRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void rampingRadiation::correct()
{}


tmp<volScalarField> rampingRadiation::Shs()
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

        if(time > rampStartTime_ + rampTimeInterval_ ){
            QrConst_.internalField() += rampStep_;
            rampStartTime_ += rampTimeInterval_;
            Info << "time, QrConst = " << time << ", " << gMax(QrConst_) << endl;
        }

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

tmp<volScalarField> rampingRadiation::ShsConst() const
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
