/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "mutUSpaldingBlowingWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutUSpaldingBlowingWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const scalarField& muw = turbModel.mu().boundaryField()[patchi];

    //const scalarField& phiw = turbModel.phi().boundaryField()[patchi];
    const scalarField& phiw = patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    scalarField fuelMassFlux( - phiw/patch().magSf()*rPhi_*1000.0); //convert to g/m2/s, and back to pyrolysate

    return max
    (
        scalar(0),
//        rhow*sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - muw
        rhow*sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL)
      * min 
        (
            scalar(1), 
            (fuelMassFlux/hOCp_+SMALL) / (Foam::exp(fuelMassFlux/scalar(hOCp_))-scalar(1) + SMALL)
        )
      - muw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutUSpaldingBlowingWallFunctionFvPatchScalarField::
mutUSpaldingBlowingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutUSpaldingWallFunctionFvPatchScalarField(p, iF),
    hOCp_(10.0),
    rPhi_(2.5)
{}


mutUSpaldingBlowingWallFunctionFvPatchScalarField::
mutUSpaldingBlowingWallFunctionFvPatchScalarField
(
    const mutUSpaldingBlowingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutUSpaldingWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    hOCp_(ptf.hOCp_),
    rPhi_(ptf.rPhi_)
{}


mutUSpaldingBlowingWallFunctionFvPatchScalarField::
mutUSpaldingBlowingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutUSpaldingWallFunctionFvPatchScalarField(p, iF, dict),
    hOCp_(dict.lookupOrDefault<scalar>("hOCp", 10.0)),
    rPhi_(dict.lookupOrDefault<scalar>("rPhi", 2.5))
{}


mutUSpaldingBlowingWallFunctionFvPatchScalarField::
mutUSpaldingBlowingWallFunctionFvPatchScalarField
(
    const mutUSpaldingBlowingWallFunctionFvPatchScalarField& wfpsf
)
:
    mutUSpaldingWallFunctionFvPatchScalarField(wfpsf),
    hOCp_(wfpsf.hOCp_),
    rPhi_(wfpsf.rPhi_)
{}


mutUSpaldingBlowingWallFunctionFvPatchScalarField::
mutUSpaldingBlowingWallFunctionFvPatchScalarField
(
    const mutUSpaldingBlowingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutUSpaldingWallFunctionFvPatchScalarField(wfpsf, iF),
    hOCp_(wfpsf.hOCp_),
    rPhi_(wfpsf.rPhi_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutUSpaldingBlowingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeKeyword("hOCp") << hOCp_ << token::END_STATEMENT << nl;
    os.writeKeyword("rPhi") << rPhi_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutUSpaldingBlowingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
