/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "greyDiffusiveRadiationFireMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "mappedPatchBase.H"
#include "mapDistribute.H"

#include "fvDOM.H"
#include "constants.H"
#include "basicSolidThermo.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

template<>
const char*
NamedEnum
<
    greyDiffusiveRadiationFireMixedFvPatchScalarField::emissivityMode,
    2
>::names[] =
{
    "lookup",
    "solidThermo"
};

const NamedEnum
<
    greyDiffusiveRadiationFireMixedFvPatchScalarField::emissivityMode,
    2
>
greyDiffusiveRadiationFireMixedFvPatchScalarField::emissivityModeNames_;


} //End namespace radiation

} //End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
greyDiffusiveRadiationFireMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefinedT"),
    emissivity_(p.size(), 0.0),
    emissivityMode_(LOOKUP)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
greyDiffusiveRadiationFireMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFireMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    emissivityMode_(ptf.emissivityMode_)
{}


Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
greyDiffusiveRadiationFireMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookup("T")),
    emissivity_(p.size(), 0.0),
    emissivityMode_(emissivityModeNames_.read(dict.lookup("emissivityMode")))
{

    if (emissivityMode_ == LOOKUP)
    {
         emissivity_ = scalarField("emissivity", dict, p.size());
    }
    else if(emissivityMode_ == SOLIDTHERMO)
    {
        if (!isA<mappedPatchBase>(this->patch().patch()))
        {
            FatalErrorIn
            (
                "greyDiffusiveRadiationFireMixedFvPatchScalarField::"
                "greyDiffusiveRadiationFireMixedFvPatchScalarField\n"
                "(\n"
                "    const fvPatch& p,\n"
                "    const DimensionedField<scalar, volMesh>& iF,\n"
                "    const dictionary& dict\n"
                ")\n"
            )   << "\n    patch type '" << p.type()
                << "' not type '" << mappedPatchBase::typeName << "'"
                << "\n    for patch " << p.name()
                << " of field " << dimensionedInternalField().name()
                << " in file " << dimensionedInternalField().objectPath()
                << exit(FatalError);
        }

        if
        (
            db().foundObject<basicSolidThermo>
            (
                "solidThermophysicalProperties"
            )
        )
        {
            FatalErrorIn
            (
                "greyDiffusiveRadiationFireMixedFvPatchScalarField::"
                "(const fvPatch& , "
                "const DimensionedField<scalar, volMesh>& "
                "const dictionary&) "
            )   << " solidThermophysicalProperties is not found "
                << " in " << db().name()
                << nl << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "greyDiffusiveRadiationFireMixedFvPatchScalarField::"
            "(const fvPatch& , "
            "const DimensionedField<scalar, volMesh>& "
            "const dictionary&) "
        )   << " mode " << emissivityMode_
            << " is not in the table " << emissivityModeNames_.toc()
            << nl << exit(FatalError);
    }

    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        const scalarField& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        refValue() =
            emissivity_*4.0*constant::physicoChemical::sigma.value()*pow4(Tp)
           /constant::mathematical::pi;

        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
greyDiffusiveRadiationFireMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFireMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    emissivityMode_(ptf.emissivityMode_)
{}


Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
greyDiffusiveRadiationFireMixedFvPatchScalarField
(
    const greyDiffusiveRadiationFireMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    emissivityMode_(ptf.emissivityMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label patchI = patch().index();

    if(emissivityMode_ == SOLIDTHERMO)
    {
        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            patch().patch()
        );
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];

        // Force recalculation of mapping and schedule
        const mappedPatchBase& distMap = mpp;

        emissivity_ =
            nbrPatch.lookupPatchField<volScalarField, scalar>("emissivity");

        mpp.distribute(emissivity_);
//        mapDistribute::distribute
//        (
//            Pstream::defaultCommsType,
//            distMap.schedule(),
//            distMap.constructSize(),
//            distMap.subMap(),           // what to send
//            distMap.constructMap(),     // what to receive
//            emissivity_
//        );
    }

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);


    if (dom.nLambda() != 1)
    {
        FatalErrorIn
        (
            "Foam::radiation::"
            "greyDiffusiveRadiationFireMixedFvPatchScalarField::updateCoeffs"
        )   << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;
    vectorField n = patch().Sf()/patch().magSf();

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    ray.Qr().boundaryField()[patchI] += Iw*(n & ray.dAve());

    forAll(Iw, faceI)
    {
        scalar Ir = 0.0;

        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& d = dom.IRay(rayI).d();

            const scalarField& IFace =
                dom.IRay(rayI).ILambda(lambdaId).boundaryField()[patchI];

            if ((-n[faceI] & d) < 0.0)
            {
                // q into the wall
                const vector& dAve = dom.IRay(rayI).dAve();
                Ir += IFace[faceI]*mag(n[faceI] & dAve);
            }
        }

        const vector& d = dom.IRay(rayId).d();

        if ((-n[faceI] & d) > 0.0)
        {
            // direction out of the wall
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] =
                (
                    Ir*(1.0 - emissivity_[faceI])
                  + emissivity_[faceI]
                  * constant::physicoChemical::sigma.value()*pow4(Tp[faceI])
                )
               /constant::mathematical::pi;

            // Emmited heat flux from this ray direction
            ray.Qem().boundaryField()[patchI][faceI] =
                refValue()[faceI]*(n[faceI] & ray.dAve());

        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used


            // Incident heat flux on this ray direction
            ray.Qin().boundaryField()[patchI][faceI] =
                Iw[faceI]*(n[faceI] & ray.dAve());
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveRadiationFireMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivityMode") << emissivityModeNames_[emissivityMode_]
        << token::END_STATEMENT << nl;
    emissivity_.writeEntry("emissivity", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveRadiationFireMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
