/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenCFD Ltd.
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

#include "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "turbulenceModel.H"
#include "thermoSingleLayer.H"
#include "pyrolysisModel.H"

#include "constants.H"

#include "radiationModel.H"
#include "absorptionEmissionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// template<>
// const char*
// NamedEnum
// <
//     turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
//     operationMode,
//     4
// >::names[] =
// {
//     "radiative_flux_from_neighbouring_region",
//     "radiative_flux_from_this_region",
//     "no_radiation_contribution",
//     "unknown"
// };

// const NamedEnum
// <
//     turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
//     operationMode,
//     4
// >
// turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
// operationModeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::filmModelType&
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }


    FatalErrorIn
    (
        "const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
        "filmModelType& "
        "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
        "filmModel() const"
    )
        << "Unable to locate film region " << filmRegionName_
        << abort(FatalError);

    return **models.begin();
}


const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
pyrolysisModelType&
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    FatalErrorIn
    (
        "const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
        "pyrolysisModelType& "
        "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
        "pyrModel() const"
    )
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << abort(FatalError);

    return **models.begin();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    KName_("undefined-K"),
    convectiveScaling_(1.0),
    convectiveCoefficient_(1.0),
    filmDeltaDry_(0.0),
    filmDeltaWet_(0.0),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
(
    const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        p,
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_),
    oldMode_(psf.oldMode_)
{}


turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    radiationCoupledBase(p, dict),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    KName_(dict.lookup("K")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    convectiveCoefficient_(dict.lookupOrDefault<scalar>("convectiveCoefficient", 1.0)),
    filmDeltaDry_
    (
            dict.lookupOrDefault<scalar>("filmDeltaDry", 0.0000)
    ),
    filmDeltaWet_
    (
            dict.lookupOrDefault<scalar>("filmDeltaWet", 0.0002)
    ),
//    filmDeltaDry_(readScalar(dict.lookup("filmDeltaDry"))),
//    filmDeltaWet_(readScalar(dict.lookup("filmDeltaWet"))),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
            "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField\n"
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


turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
(
    const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField&
        psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        psf.patch(),
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    KName_(psf.KName_),
    convectiveScaling_(psf.convectiveScaling_),
    convectiveCoefficient_(psf.convectiveCoefficient_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_),
    oldMode_(psf.oldMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();
    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    scalarField intFld(patchInternalField());

    const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField&
        nbrField =
        refCast
        <
            const turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField& Tp = *this;

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    scalarField nbrConvFlux(KDeltaNbr*(intFld - nbrIntFld));

    scalarField nbrTotalFlux = nbrConvFlux;
    // scalarList nbrRadField(nbrPatch.size(), 0.0);
    scalarList QrCoupled(nbrPatch.size(), 0.0);//kvm
    scalarList Tfilm(nbrPatch.size(), 0.0);//kvm
    scalarList alpha(nbrPatch.size(), 0.0);//kvm
    scalarList filmConv(nbrPatch.size(), 0.0);//kvm
    scalarList filmDelta(nbrPatch.size(), 0.0);//kvm
    scalarList myRadField(patch().size(), 0.0);

    const pyrolysisModelType& pyrolysis = pyrModel();
    const filmModelType& film = filmModel();

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {

        scalarField nbrConvFlux(convectiveCoefficient_*(intFld - nbrIntFld));
        
        const label filmPatchI =
            pyrolysis.nbrCoupledPatchID(film, patchI);

        const scalarField Qconvw(film.Qconvw(filmPatchI));

        // kvm, Qconvw is not right
        filmConv =
            // pyrRegion.mapRegionPatchField
            // (
            //     filmModel,
            //     patchI,
            //     filmPatchI,
            //     Qconvw,
            //     true
            // );
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "qFilmToWall",
                patchI,
                true
            );

        QrCoupled =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Qin",
                patchI,
                true
            );

        // Info << "QrCoupled " << QrCoupled << endl;

        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                patchI,
                true
            );

        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                patchI,
                true
            );
        alpha =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "alpha",
                patchI,
                true
            );

            const fvMesh& mesh = patch().boundaryMesh().mesh();
            if (! (mesh.foundObject<radiation::radiationModel>("radiationProperties")))
            {
                FatalErrorIn
                (
                    "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::"
                    "turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const DimensionedField<scalar, volMesh>& iF,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n    radiationProperties file not found in pyrolysis region\n" 
                    << exit(FatalError);
            }
            const radiation::radiationModel& radiation =
                mesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField temissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );
            

        // nbrRadField =
        //     nbrPatch.lookupPatchField<volScalarField, scalar>
        //     (
        //         neighbourFieldRadiativeName_
        //     );

        // Note: the Qr radiative flux is positive outgoing.
        // For a hot solid radiating into a cold fluid Qr will be negative.


        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        // mpp.distribute(nbrRadField);

    //Qrad is negative going out of the solid
    //Qcon is positive going out of the solid

        // emissivity_ =
        //     patch().lookupPatchField<volScalarField, scalar>("emissivity");


        if (debug)
        {
            // scalar Qr = gSum(nbrRadField*patch().magSf());

            Info<< mesh.name() << ':'
                << patch().name() << ':'
                << this->dimensionedInternalField().name() << " :" << nl
                // << "     radiative heat  [W] : " << Qr << nl
//kvm                << "     predicted wallT [K] : " << gAverage(Twall) << nl
                << endl;
        }

        label nFixed = 0;

         // Estimate wetness of the film (1: wet , 0: dry)
         scalarField ratio
         (
            min
            (
                max
                (
                    (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                    scalar(0.0)
                ),
                scalar(1.0)
            )
         );

        forAll(*this, i)
        {
            scalar qConvWet = -filmConv[i];
            scalar qConvDry = nbrConvFlux[i];

            const scalar sigma = constant::physicoChemical::sigma.value();
            scalar qRadWet = 0.0; // all film absorption takes place in film model
            // use of 'ratio' leads to non-conservation of energy
            // if use 'ratio', then need to use it in film model absorption as well
            //kvm scalar qRadWet =-QrCoupled[i];
            scalar qRadDry =
               -temissivity[i]*QrCoupled[i]
               //kvm -(1.0-ratio)*temissivity[i]*QrCoupled[i]
               +temissivity[i]*sigma*pow4(operator[](i));

            //DEBUG(temissivity[i]);
            //DEBUG(QrCoupled[i]);
            //DEBUG(operator[](i));

            scalar qConv = alpha[i]*qConvWet + (1.0-alpha[i])*qConvDry;
            scalar qRad  = alpha[i]*qRadWet  + (1.0-alpha[i])*qRadDry;

            //DEBUG(qConvDry);
            //DEBUG(qConvWet);
            //DEBUG(qRadDry);
            //DEBUG(qRadWet);
            //prateep's scaling for 16 cells across flue space is 1.0
            //                       8 cells across flue space is 2.0
            //                       4 cells across flue space is 4.0 
            qConv *= convectiveScaling_;

            nbrTotalFlux[i] = qConv + qRad;
            this->refValue()[i] = operator[](i);  // not used
            this->refGrad()[i] = -nbrTotalFlux[i]/kappa(*this)()[i];                
            //DEBUG(nbrTotalFlux[i]);
            //DEBUG(kappa(*this)()[i]);

            this->valueFraction()[i] = 0.0;
            nFixed++;
        }

        if (debug)
        {
            Pout<< "Using " << nFixed << " fixedValue out of " << this->size()
                << endl;
        }
    }
    else // In fluid
    {
//        typedef regionModels::pyrolysisModels::pyrolysisModel
//            pyrolysisModelType;
//
//        const regionModels::regionModel& pyrolysisRegion =
//            db().time().lookupObject<regionModels::regionModel>
//            (
//                "pyrolysisProperties"
//            );
//        const pyrolysisModelType& pyrolysisModel =
//            dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
//
//        pyrolysisModelType& pyrolysis =
//            const_cast<pyrolysisModelType&>(pyrolysisModel);
//
//        const regionModels::regionModel& filmModel =
//            pyrolysisModel.db().time().lookupObject<regionModels::regionModel>
//            (
//                "surfaceFilmProperties"
//            );

        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                nbrPatchI,
                true
            );

        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                nbrPatchI,
                true
            );

        const scalarField& myRadField =
            patch().lookupPatchField<volScalarField, scalar>
            (
                fieldRadiativeName_
            );

        // do we still need to do mpp.distribute Tfilm
        mpp.distribute(Tfilm);

        mpp.distribute(filmDelta);

        // use solid internal cell as Twall of gas. See if we can make re-radiation stable.
        // kvm, superseded by below  scalarField Twall = nbrIntFld;
        scalarList Twall(patch().size(), 0.0);//kvm
        
        // Estimate wetness of the film (1: wet , 0: dry)
        scalarField ratio
        (
           min
           (
               max
               (
                   (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                   scalar(0.0)
               ),
               scalar(1.0)
           )
        );

        forAll(*this, i)
        {
            scalar filmDeltaWet=0.0002; //m
            scalar filmDeltaDry=0.0000; //m

            scalar Twet = min(max(Tfilm[i], 298.15), 378.4);
            scalar Tdry = nbrIntFld[i];

            Twall[i] = ratio[i]*(Twet - Tdry) + Tdry;
        }

        if (debug)
        {
            scalar Qr = gSum(myRadField*patch().magSf());

            // Info<< mesh.name() << ':'
            //     << patch().name() << ':'
            //     << this->dimensionedInternalField().name() << " :" << nl
            //     << "     radiative heat  [W] : " << Qr << nl
            //     << "     predicted wallT [K] : " << gAverage(Twall) << nl
            //     << endl;
        }

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qc = gSum(nbrConvFlux*patch().magSf());
        // scalar Qr = gSum(radField*patch().magSf());
        scalar Qt = gSum(nbrTotalFlux*patch().magSf());

        // Info<< mesh.name() << ':'
        //     << patch().name() << ':'
        //     << this->dimensionedInternalField().name() << " -> "
        //     << nbrMesh.name() << ':'
        //     << nbrPatch.name() << ':'
        //     << this->dimensionedInternalField().name() << " :"
        //     << " heatFlux:" << Qc
        //     // << " radiativeFlux:" << Qr
        //     << " totalFlux:" << Qt
        //     << " walltemperature "
        //     << " min:" << gMin(*this)
        //     << " max:" << gMax(*this)
        //     << " avg:" << gAverage(*this)
        //     << endl;
    }
}


void turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os,
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    writeEntryIfDifferent<word>
    (
        os,
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("K")<< KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<<
        fieldRadiativeName_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
    radiationCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
