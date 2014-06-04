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

#include "constHTemperatureRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::constHTemperatureRadiationFvPatchScalarField::
constHTemperatureRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    Tinf_(p.size(), 0.0),
    emissivity_(p.size(),0.0),
    h_(p.size(), 0.0),
    modeName_("undefined-modeName")
{
    refValue() = 295;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::constHTemperatureRadiationFvPatchScalarField::
constHTemperatureRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    Tinf_("Tinf", dict, p.size()),
    emissivity_("emissivity", dict, p.size()),
    h_("h", dict, p.size()),
    modeName_(dict.lookup("mode"))
{
    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}

Foam::constHTemperatureRadiationFvPatchScalarField::
constHTemperatureRadiationFvPatchScalarField
(
    const constHTemperatureRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    Tinf_(ptf.Tinf_, mapper),
    emissivity_(ptf.emissivity_, mapper),
    h_(ptf.h_, mapper),
    modeName_(ptf.modeName_)
{}


Foam::constHTemperatureRadiationFvPatchScalarField::
constHTemperatureRadiationFvPatchScalarField
(
    const constHTemperatureRadiationFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    Tinf_(tppsf.Tinf_),
    emissivity_(tppsf.emissivity_),
    h_(tppsf.h_),
    modeName_(tppsf.modeName_)
{}

Foam::constHTemperatureRadiationFvPatchScalarField::
constHTemperatureRadiationFvPatchScalarField
(
    const constHTemperatureRadiationFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    Tinf_(tppsf.Tinf_),
    emissivity_(tppsf.emissivity_),
    h_(tppsf.h_),
    modeName_(tppsf.modeName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constHTemperatureRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Tinf_.autoMap(m);
    emissivity_.autoMap(m);
    h_.autoMap(m);
}


void Foam::constHTemperatureRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const constHTemperatureRadiationFvPatchScalarField& tiptf =
         refCast<const constHTemperatureRadiationFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
    emissivity_.rmap(tiptf.emissivity_, addr);
    h_.rmap(tiptf.h_, addr);
}

void Foam::constHTemperatureRadiationFvPatchScalarField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    // const scalarField K_ = patch().lookupPatchField<volScalarField, scalar>("K");
    
    const scalar sigma = constant::physicoChemical::sigma.value();
    
    if( modeName_ == "fixed" ){
        // do nothing
    }
    else if( modeName_ == "correlation" ){
        // TODO: get all these properties from the gas phase
        scalar g=9.8;
        scalar Tw=operator[](0);
        scalarField& Tfaces = *this;
        scalar Tamb=average(Tfaces);
        scalar Tf=0.5*(Tw+Tamb);
        scalar beta=1.0/Tf;
        scalar L=0.9;
        scalar alpha=38.3e-6;
        scalar nu=26.4e-6;
        
        scalar RaL = g*beta*(Tw-Tamb)*pow(L,3)/(nu*alpha);
        // Negative RaL cause corellation to blow up numerically
        RaL = max(RaL,0.0);

        scalar Pr=0.69;
        scalar NuL=pow(
            0.825+(0.387)*pow(RaL,1./6.)
             /
            pow(
                1.0+pow(
                    0.492/Pr
                    ,9./16.)
                ,8./27.)
            ,2.);
        
        scalar k=33.8e-3;
        dimensionedScalar hbar("hbar",dimEnergy/dimArea/dimTime/dimTemperature,NuL*k/L);
        // TODO: something is wrong with this bc.  Tamb probably not right.
        // Info << "hTamb " << Tamb << endl;
        // Info << "hTw " << Tw << endl;
        // Info << "hbeta " << beta << endl;
        // Info << "hbetaDT " << beta*(Tw-Tamb) << endl;
        // Info << "hTf " << Tf << endl;
        // Info << "RaL " << RaL << endl;
        // Info << "NuL " << NuL << endl;
        // Info << "hbar " << hbar << endl;
        forAll(h_,i){
            h_[i]=hbar.value();
        } 
    }
    else
    {
        FatalErrorIn
        (
            "constHTemperatureRadiationFvPatchScalarField"
            "::updateCoeffs()"
            )   << "Illegal mode, " 
                << "choose from:\n"
                << "\tfixed\n"
                << "\tcorrelation"
                << exit(FatalError);
    }

    forAll(*this, i)
    {
        const scalar T = operator[](i);

        // positive heat flux heats solid, negative cools solid
        // convection
        const scalar qConv = h_[i]*(Tinf_[i] - T);

        // radiation
        const scalar absorptivity = emissivity_[i];
        const scalar emissivitySurroundings = emissivity_[i];
        const scalar qRadIncident = emissivitySurroundings*sigma*pow4(Tinf_[i]);
        const scalar qRadAbsorbed = + absorptivity  *qRadIncident;
        const scalar qRadEmission = - emissivity_[i]*sigma*pow4(T);
        const scalar qRad = 
            + qRadAbsorbed 
            + qRadEmission;
        const scalar qTotal = qRad + qConv;

        this->refValue()[i] = Tinf_[i];  // not used
        this->refGrad()[i] = qTotal;
        this->valueFraction()[i] = 0.0;
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::constHTemperatureRadiationFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    Tinf_.writeEntry("Tinf", os);
    emissivity_.writeEntry("emissivity", os);
    h_.writeEntry("h", os);
    os.writeKeyword("mode")<< modeName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        constHTemperatureRadiationFvPatchScalarField
    );

}

// ************************************************************************* //
