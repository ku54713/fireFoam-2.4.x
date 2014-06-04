/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "LESModel.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::
turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    KName_("undefined-K"),
    emissivity_(p.size(), 0.0),
    cQcScale_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::
turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourFieldRadiativeName_(ptf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(ptf.fieldRadiativeName_),
    KName_(ptf.KName_),
    emissivity_(ptf.emissivity_),
    cQcScale_(ptf.cQcScale_)
{}


turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::
turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    KName_(dict.lookup("K")),
    emissivity_(p.size(), 0.0),
    cQcScale_(readScalar(dict.lookup("cQcScale")))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::"
            "turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField\n"
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

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::
turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField&
        wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourFieldRadiativeName_(wtcsf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(wtcsf.fieldRadiativeName_),
    KName_(wtcsf.KName_),
    emissivity_(wtcsf.emissivity_),
    cQcScale_(wtcsf.cQcScale_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KName_ == "none")
    {
        const compressible::LESModel& model =
            db().lookupObject<compressible::LESModel>("LESProperties");

        tmp<volScalarField> talpha = model.alphaEff();

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            talpha().boundaryField()[patch().index()]
           *thermo.Cp()().boundaryField()[patch().index()];
    }
    else if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
    {
        return patch().lookupPatchField<volScalarField, scalar>(KName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KName_))
    {
        const symmTensorField& KWall =
            patch().lookupPatchField<volSymmTensorField, scalar>(KName_);

        vectorField n = patch().nf();

        return n & KWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField"
            "::K() const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

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
    //const mapDistribute& distMap = mpp.map();

    scalarField intFld = patchInternalField();

    const turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField&
    nbrField =  refCast
        <const turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld = nbrField.patchInternalField();
    mpp.distribute(nbrIntFld);
    //mapDistribute::distribute
    //(
    //    Pstream::defaultCommsType,
    //    distMap.schedule(),
    //    distMap.constructSize(),
    //    distMap.subMap(),           // what to send
    //    distMap.constructMap(),     // what to receive
    //    nbrIntFld
    //);

    // Swap to obtain full local values of neighbour K*delta
    scalarField nbrKDelta = nbrField.K()*nbrPatch.deltaCoeffs();
    mpp.distribute(nbrKDelta);
    //mapDistribute::distribute
    //(
    //    Pstream::defaultCommsType,
    //    distMap.schedule(),
    //    distMap.constructSize(),
    //    distMap.subMap(),           // what to send
    //    distMap.constructMap(),     // what to receive
    //    nbrKDelta
    //);

    //scalarField myKDelta = K()*patch().deltaCoeffs();
    //myKDelta = K();

    //scalarField nbrConvFlux = nbrKDelta*(*this - nbrIntFld);
    scalarField nbrConvFlux = cQcScale_*nbrKDelta*(intFld - nbrIntFld);

    scalarField nbrTotalFlux = nbrConvFlux;

    scalarList radField(nbrPatch.size(),0.0);  
    scalarField Twall(patch().size(),0.0);

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {
        radField =
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldRadiativeName_
            );

        // Swap to obtain full local values of neighbour radiative heat flux field
        mpp.distribute(radField);
        //mapDistribute::distribute
        //(
        //    Pstream::defaultCommsType,
        //    distMap.schedule(),
        //    distMap.constructSize(),
        //    distMap.subMap(),           // what to send
        //    distMap.constructMap(),     // what to receive
        //    radField
        //);

    //Qrad is negative going out of the solid
    //Qcon is positive going out of the solid

        emissivity_ =
            patch().lookupPatchField<volScalarField, scalar>("emissivity");

        nbrTotalFlux -= radField*emissivity_;
        nbrTotalFlux += emissivity_*constant::physicoChemical::sigma.value()*pow(*this,4);

        forAll(*this, i)
        {
                //fixed gradient BC, use internal T to replace Tb for re-radiation
                this->refValue()[i] = operator[](i);  // not used
                this->refGrad()[i] = -nbrTotalFlux[i]/K()()[i];
                this->valueFraction()[i] = 0.0;
        }

    }
    else // In fluid
    {
        //bug fix: radField need to be resize to local patch size. will cause trouble in debug writing
        //const scalarField& radField =
        radField =       
            patch().lookupPatchField<volScalarField, scalar>
            (
                fieldRadiativeName_
            );

        Twall = nbrIntFld;

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qc = gSum(nbrConvFlux*patch().magSf());
        scalar Qr = gSum(radField*patch().magSf());
        scalar Qt = gSum(nbrTotalFlux*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heatFlux:" << Qc
            << " radiativeFlux:" << Qr
            << " totalFlux:" << Qt
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< fieldRadiativeName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("cQcScale") << cQcScale_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadiationQinQcScaleCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
