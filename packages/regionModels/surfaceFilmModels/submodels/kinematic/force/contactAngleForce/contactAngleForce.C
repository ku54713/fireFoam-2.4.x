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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "fvPatchField.H"
#include "patchDist.H"
#include "zeroGradientFvPatchFields.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);
addToRunTimeSelectionTable(force, contactAngleForce, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void contactAngleForce::initialise()
{
    const wordReList zeroForcePatches(coeffDict_.lookup("zeroForcePatches"));

    if (zeroForcePatches.size())
    {
        const polyBoundaryMesh& pbm = owner_.regionMesh().boundaryMesh();
        scalar dLim = readScalar(coeffDict_.lookup("zeroForceDistance"));

        Info<< "        Assigning zero contact force within " << dLim
            << " of patches:" << endl;

        labelHashSet patchIDs = pbm.patchSet(zeroForcePatches);

        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            label patchI = iter.key();
            Info<< "            " << pbm[patchI].name() << endl;
        }

        patchDist dist(owner_.regionMesh(), patchIDs);

        mask_ = pos(dist - dimensionedScalar("dLim", dimLength, dLim));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    force(typeName, owner, dict),
    Ccf_(readScalar(coeffDict_.lookup("Ccf"))),
    rndGen_(label(0), -1),
    distribution_
    (
        distributionModels::distributionModel::New
        (
            coeffDict_.subDict("contactAngleDistribution"),
            rndGen_
        )
    ),
    timeIntervalDistribution_
    (
        distributionModels::distributionModel::New
        (
            coeffDict_.subDict("timeIntervalDistribution"),
            rndGen_
        )
    ),
    mask_
    (
        IOobject
        (
            typeName + "_contactForceMask",
            owner_.time().timeName(),
            owner_.regionMesh(),
            // IOobject::MUST_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        owner_.regionMesh(),
        dimensionedScalar("mask", dimless, 1.0)
        // zeroGradientFvPatchScalarField::typeName
     ),
    contactAngle_
    (
        IOobject
        (
            typeName + "_contactAngle",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    contactAngleOld_
    (
        IOobject
        (
            typeName + "_contactAngleOld",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    contactAngleNew_
    (
        IOobject
        (
            typeName + "_contactAngleNew",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    timeOld_
    (
        IOobject
        (
            typeName + "_timeOld",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    timeInterval_
    (
        IOobject
        (
            typeName + "_timeInterval",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
     ),
    nHits_
    (
        IOobject
        (
            typeName + "_nHits",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
     )
{
    initialise();
    // mask for zero-ing out contact angle
    // mask_ = pos(mask_-0.5);
    forAll(contactAngleOld_,cellI){
        contactAngleOld_[cellI] = distribution_->sample();
        timeOld_[cellI] = owner.time().value();
    }
    forAll(contactAngleNew_,cellI){
        contactAngleNew_[cellI] = distribution_->sample();
        timeInterval_[cellI] = timeIntervalDistribution_->sample();
    }
    const kinematicSingleLayer& film =
        dynamic_cast<const kinematicSingleLayer&>(owner_);
    Info << "dx " << sqrt(film.magSf()[0]) << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

contactAngleForce::~contactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> contactAngleForce::correct(volVectorField& U)
{
    const kinematicSingleLayer& film =
        dynamic_cast<const kinematicSingleLayer&>(owner_);

    tmp<volVectorField> tForce
    (
        new volVectorField
        (
            IOobject
            (
                typeName + "_contactForce",
                owner_.time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner_.regionMesh(),
            dimensionedVector("zero", dimForce/dimArea, vector::zero),
            // eliminate contact-line at boundaries
            zeroGradientFvPatchVectorField::typeName
        )
    );

    scalar time=film.time().value();

    forAll(contactAngle_,cellI){
        if(time > timeOld_[cellI]+timeInterval_[cellI]){
            contactAngleOld_[cellI]=contactAngleNew_[cellI];
            timeOld_[cellI]=timeOld_[cellI]+timeInterval_[cellI];
            contactAngleNew_[cellI]=distribution_->sample();
            timeInterval_[cellI]=timeIntervalDistribution_->sample();
        }
        scalar f = (time - timeOld_[cellI])/(timeInterval_[cellI]);
        contactAngle_[cellI] = (1.0-f)*contactAngleOld_[cellI]+f*contactAngleNew_[cellI];
    }


    contactAngleNew_.correctBoundaryConditions();
    contactAngleOld_.correctBoundaryConditions();
    contactAngle_.correctBoundaryConditions();
    timeInterval_.correctBoundaryConditions();

    vectorField& force = tForce().internalField();

    const labelUList& own = owner_.regionMesh().owner();
    const labelUList& nbr = owner_.regionMesh().neighbour();

    const scalarField& magSf = owner_.magSf();

    const volScalarField& alpha = owner_.alpha();
    const volScalarField& sigma = owner_.sigma();
    const volScalarField& deltaf = film.delta();

    volVectorField gradAlpha(fvc::grad(alpha));

    scalarField nHits(film.regionMesh().nCells(), 0.0);

    // minimum film thickness used for scaling contact angle force
    const scalar deltaf0 = 2e-4;

    // list of all faces, with own and nbr being the cell pair
    forAll(nbr, faceI)
    {
        const label cellO = own[faceI];
        const label cellN = nbr[faceI];

        label cellI = -1;
        label cellOther = -1;
        if ((alpha[cellO] > 0.5) && (alpha[cellN] < 0.5))
        {
            cellI = cellO;
            cellOther = cellN;
        }
        else if ((alpha[cellO] < 0.5) && (alpha[cellN] > 0.5))
        {
            cellI = cellN;
            cellOther = cellO;
        }

        // const scalar dxInverse = film.regionMesh().deltaCoeffs()[faceI];
        // const scalar area = film.regionMesh().magSf()[faceI];
        if (cellI != -1)
        {
            // is this the right length to be using?
            const scalar dxInverse = film.regionMesh().deltaCoeffs()[faceI]; // really 1/dxInverse
            // const scalar area = film.regionMesh().magSf()[faceI]; // really 1/dxInverse
            const vector n =
                gradAlpha[cellI]/(mag(gradAlpha[cellI]) + ROOTVSMALL);
            // won't the contact angle be changing each timestep?
            // scalar theta = cos(degToRad(distribution_->sample()));
            scalar ratio = min(deltaf[cellI]/deltaf0,1.0);
            scalar theta = cos(contactAngle_[cellI]);
            // force[cellI] += mask_[cellI]*Ccf_*n*sigma[cellI]*(1.0 - theta)/dxInverse;
            force[cellI] += mask_[cellI]*n*sigma[cellI]*(1.0 - theta)/Ccf_*ratio; // using Ccf_ as characteristic length
            // TODO: seems like we want to have the contact angel force on both sides of
            // the contact line interface.  This way, the residual liquid will be forced out.
            // (not sure about this one)
            // force[cellOther] += mask_[cellI]*Ccf_*n*sigma[cellI]*(1.0 - theta)/dxInverse;
            // force[cellI] += mask_[cellI]*Ccf_*n*sigma[cellI]*(1.0 - theta)*dxInverse; //bug fix
            // force[cellI] += mask_[cellI]*Ccf_*n*sigma[cellI]*(1.0 - theta)*area; //bug fix
            nHits[cellI]++;
        }
    }

    forAll(alpha.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& alphaf = alpha.boundaryField()[patchI];
        const scalarField& dxInverse = alphaf.patch().deltaCoeffs();
        // const scalarField& area = alphaf.patch().magSf();
        // why do we loop over the top and bottom patches?
        const labelUList& faceCells = alphaf.patch().faceCells();

        forAll(alphaf, faceI)
        {
            label cellO = faceCells[faceI];

            if ((alpha[cellO] > 0.5) && (alphaf[faceI] < 0.5))
            {
                // won't the contact angle be changing each timestep?
                const vector n =
                    gradAlpha[cellO]/(mag(gradAlpha[cellO]) + ROOTVSMALL);
                // scalar theta = cos(degToRad(distribution_->sample()));
                // ratio is important for code stability
                // when alpha=1 and deltaf is very small, the contact angle force
                // can be so large as to cause instabilities in the velocity.
                scalar ratio = min(deltaf[cellO]/deltaf0,1.0);
                scalar theta = cos(contactAngle_[cellO]);
                // is this the right dxInverse to be using?
                // force[cellO] += mask_[cellO]*Ccf_*n*sigma[cellO]*(1.0 - theta)/dxInverse[faceI];
                force[cellO] += mask_[cellO]*n*sigma[cellO]*(1.0 - theta)/Ccf_*ratio;
                // force[cellO] += mask_[cellO]*Ccf_*n*sigma[cellO]*(1.0 - theta)*area[faceI];
                nHits[cellO]++;
            }
        }
    }

    // why do we divide by nHits?  doesn't make sense
    // force /= (max(nHits, scalar(1.0))*magSf);
    // force /= magSf;
    tForce().correctBoundaryConditions();

    if (owner_.regionMesh().time().outputTime())
    {
        tForce().write();
        nHits_.internalField() = nHits;
    }

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume)); // why is this *dimVolume?

    tfvm() += tForce; // doesn't this need to be *volume?

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
