/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "LESModel.H"

#include "constants.H"

#define DEBUG(x) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

template<>
const char*
NamedEnum
<
    turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
    operationMode,
    4
>::names[] =
{
    "radiative_flux_from_neighbouring_region",
    "radiative_flux_from_this_region",
    "no_radiation_contribution",
    "unknown"
};

const NamedEnum
<
    turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
    operationMode,
    4
>
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
operationModeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    neighbourFieldConvectiveName_("undefined-neigbourFieldConvectiveName"),
    fieldConvectiveName_("undefined-fieldConvectiveName"),
    KName_("undefined-K"),
    convectiveScaling_(1.0),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
(
    const turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourFieldRadiativeName_(ptf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(ptf.fieldRadiativeName_),
    neighbourFieldConvectiveName_(ptf.neighbourFieldConvectiveName_),
    fieldConvectiveName_(ptf.fieldConvectiveName_),
    KName_(ptf.KName_),
    convectiveScaling_(ptf.convectiveScaling_),
    emissivity_(ptf.emissivity_),
    oldMode_(ptf.oldMode_)
{}


turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
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
    neighbourFieldConvectiveName_(dict.lookup("neighbourFieldConvectiveName")),
    fieldConvectiveName_(dict.lookup("fieldConvectiveName")),
    KName_(dict.lookup("K")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    emissivity_(p.size(), 0.0),
    oldMode_(unknown)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::"
            "turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField\n"
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
    Info << "convective scaling set to " << convectiveScaling_ << endl;
}


turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
(
    const turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField&
        wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourFieldRadiativeName_(wtcsf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(wtcsf.fieldRadiativeName_),
    neighbourFieldConvectiveName_(wtcsf.neighbourFieldConvectiveName_),
    fieldConvectiveName_(wtcsf.fieldConvectiveName_),
    KName_(wtcsf.KName_),
    convectiveScaling_(wtcsf.convectiveScaling_),
    emissivity_(wtcsf.emissivity_),
    oldMode_(wtcsf.oldMode_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KName_ == "none")
    {
        const compressible::LESModel& model =
            db().lookupObject<compressible::LESModel>("LESProperties");

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            model.alphaEff()().boundaryField()[patch().index()]
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
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::K()"
            " const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::
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

    scalarField intFld = patchInternalField();

    const turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField&
    nbrField =
        refCast
        <
        const turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld = nbrField.patchInternalField();
    mpp.distribute(nbrIntFld);
//    mapDistribute::distribute
//    (
//        Pstream::defaultCommsType,
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrIntFld
//    );

    if (debug)
    {
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << "  internalT "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;

        Info<< nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << "  internalT "
            << " min:" << gMin(nbrIntFld)
            << " max:" << gMax(nbrIntFld)
            << " avg:" << gAverage(nbrIntFld)
            << endl;
    }



    // Check how to operate
    operationMode mode = unknown;
    {
        if (neighbourFieldRadiativeName_ != "none")
        {
            if
            (
                nbrMesh.foundObject<volScalarField>
                (
                    neighbourFieldRadiativeName_
                )
            )
            {
                mode = radFromNeighbour;
            }
            else
            {
                mode = noRad;
            }
        }
        else
        {
            if
            (
                patch().boundaryMesh().mesh().foundObject<volScalarField>
                (
                    fieldRadiativeName_
                )
            )
            {
                mode = radFromMe;
            }
            else
            {
                mode = noRad;
            }
        }

        // Do some warnings if change of mode.
        if (0&&mode != oldMode_)
        {
            WarningIn
            (
                "turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField"
                "::updateCoeffs()"
            )   << "Switched from mode " << operationModeNames[oldMode_]
                << " to mode " << operationModeNames[mode]
                << endl;
        }
        oldMode_ = mode;
    }



    // Swap to obtain full local values of neighbour K*delta
    scalarField nbrKDelta = nbrField.K()*nbrPatch.deltaCoeffs();
    mpp.distribute(nbrKDelta);
//    mapDistribute::distribute
//    (
//        Pstream::defaultCommsType,
//        distMap.schedule(),
//        distMap.constructSize(),
//        distMap.subMap(),           // what to send
//        distMap.constructMap(),     // what to receive
//        nbrKDelta
//    );

    scalarField myKDelta = K()*patch().deltaCoeffs();

    //kvm scalarField nbrConvFlux = nbrKDelta*(*this - nbrIntFld);
    scalarField nbrConvFlux = nbrKDelta*(intFld - nbrIntFld);

    scalarField nbrTotalFlux = nbrConvFlux;
    //scalarField nbrTotalFlux;    
    
    scalarList nbrRadField(nbrPatch.size(), 0.0);
    scalarList convField(nbrPatch.size(),0.0);  //YW    
    
    scalarList QrCoupled(nbrPatch.size(), 0.0);//kvm
    //const fvPatchField<scalar>& convectiveFlux(nbrPatch.size(), 0.0);//kvm
    scalarList Tfilm(nbrPatch.size(), 0.0);//kvm
    scalarList filmConv(nbrPatch.size(), 0.0);//kvm
    scalarList filmDelta(nbrPatch.size(), 0.0);//kvm
    scalarList myRadField(patch().size(), 0.0);

    // solid
    if (mode == radFromNeighbour)
    {
        /*QrCoupled is mapped from the film region and is equal to the radiation from the primary region less the amount absorbed by the film*/
        QrCoupled =
            patch().lookupPatchField<volScalarField, scalar>("Qr");//kvm
        scalarField& convectiveFlux = const_cast< fvPatchField<scalar>& >(patch().lookupPatchField<volScalarField, scalar>("convectiveFlux"));//kvm
        //fvPatchField<scalar>& convectiveFlux = const_cast< fvPatchField<scalar>& >(patch().lookupPatchField<volScalarField, scalar>("convectiveFlux"));//kvm
        Tfilm =
            patch().lookupPatchField<volScalarField, scalar>("filmTemperature");//kvm
        filmConv =
            patch().lookupPatchField<volScalarField, scalar>("filmConv");//kvm
        filmDelta =
            patch().lookupPatchField<volScalarField, scalar>("filmDelta");//kvm

        nbrRadField =

            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldRadiativeName_
            );

        // Note: the Qr radiative flux is positive outgoing.
        // For a hot solid radiating into a cold fluid Qr will be negative.


        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        mpp.distribute(nbrRadField);


        convField =
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldConvectiveName_
            );

        // Swap to obtain full local values of neighbour radiative heat flux field
        mpp.distribute(convField);        
                
        
//        mapDistribute::distribute
//        (
//            Pstream::defaultCommsType,
//            distMap.schedule(),
//            distMap.constructSize(),
//            distMap.subMap(),           // what to send
//            distMap.constructMap(),     // what to receive
//            nbrRadField
//        );

    //Qrad is negative going out of the solid
    //Qcon is positive going out of the solid

        emissivity_ =
            patch().lookupPatchField<volScalarField, scalar>("emissivity");

        //nbrTotalFlux -= nbrRadField;//kvm *emissivity_;
        //kvm,  accounted for below  nbrTotalFlux -= QrCoupled;//kvm 

        /*const scalarField Twall =*/
            /*(nbrRadField + myKDelta*intFld + nbrKDelta*nbrIntFld)*/
            /*/(myKDelta + nbrKDelta);*/


        if (debug)
        {
            scalar Qr = gSum(nbrRadField*patch().magSf());

            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->dimensionedInternalField().name() << " :" << nl
                << "     radiative heat  [W] : " << Qr << nl
//kvm                << "     predicted wallT [K] : " << gAverage(Twall) << nl
                << endl;
        }

        label nFixed = 0;

        forAll(*this, i)
        {

            if(0&&filmDelta[i]>1e-6){ //wet
                nbrTotalFlux[i] = -filmConv[i];//kvm, "negative of" because originally computed from film's perspective
                nbrTotalFlux[i] -= QrCoupled[i];//kvm, this should be zero
                /* store convective flux for diagnostic purposes */
                // flow into pyrolysis region is negative
                convectiveFlux[i] = -filmConv[i];

            }
            else if (0){ //dry
                //nbrTotalFlux[i] = nbrConvFlux[i];
                nbrTotalFlux[i] = -convField[i]; //YW
                        
                //pqr nbrTotalFlux[i] -= QrCoupled[i];//absorption
                nbrTotalFlux[i] -= emissivity_[i]*QrCoupled[i];//absorption
                nbrTotalFlux[i] += emissivity_[i]*constant::physicoChemical::sigma.value()*pow(operator[](i),4);//emission
                /* store convective flux for diagnostic purposes */
                // flow into pyrolysis region is negative
                //convectiveFlux[i] = nbrConvFlux[i];
                convectiveFlux[i] = -convField[i];    //YW            
            }
            else{
                scalar filmDeltaWet=0.0002; //m
                scalar filmDeltaDry=0.0000; //m

                scalar ratio=(filmDelta[i]-filmDeltaDry)/(filmDeltaWet-filmDeltaDry);
                ratio=max(ratio,0.0);
                ratio=min(ratio,1.0);

                scalar qConvWet=-filmConv[i];
                //scalar qConvDry=nbrConvFlux[i];
                scalar qConvDry=-convField[i];  //YW              

                scalar qRadWet =-QrCoupled[i];
                scalar qRadDry =-emissivity_[i]*QrCoupled[i]+emissivity_[i]*constant::physicoChemical::sigma.value()*pow(operator[](i),4);
                
                scalar qConv=ratio*(qConvWet-qConvDry) + qConvDry;
                scalar qRad =ratio*(qRadWet-qRadDry)   + qRadDry;

                //prateep's scaling for 16 cells across flue space is 1.0
                //                       8 cells across flue space is 2.0
                //                       4 cells across flue space is 4.0 
                //qConv*=2.0;
                qConv*=convectiveScaling_;

                convectiveFlux[i] = qConv;

                nbrTotalFlux[i] = qConv+qRad;

            }
            this->refValue()[i] = operator[](i);  // not used
            //kvm this->refValue()[i] = Twall[i];
            this->refGrad()[i] = -nbrTotalFlux[i]/K()()[i];
            //kvm this->refGrad()[i] = 0.0;   // not used
            this->valueFraction()[i] = 0.0;
            //kvm this->valueFraction()[i] = 1.0;
            nFixed++;
         }

        if (debug)
        {
            Pout<< "Using " << nFixed << " fixedValue out of " << this->size()
                << endl;
        }
    }
    else if (mode == radFromMe) //fluid
    {
        Tfilm =
            nbrPatch.lookupPatchField<volScalarField, scalar>("filmTemperature");//kvm
        filmDelta =
            nbrPatch.lookupPatchField<volScalarField, scalar>("filmDelta");//kvm

        const scalarField& myRadField =
            patch().lookupPatchField<volScalarField, scalar>
            (
             fieldRadiativeName_
            );
        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        mpp.distribute(Tfilm);
//        mapDistribute::distribute
//        (
//            Pstream::defaultCommsType,
//            distMap.schedule(),
//            distMap.constructSize(),
//            distMap.subMap(),           // what to send
//            distMap.constructMap(),     // what to receive
//            Tfilm
//        );
        // Swap to obtain full local values of neighbour radiative heat flux
        // field
        mpp.distribute(filmDelta);
//        mapDistribute::distribute
//        (
//            Pstream::defaultCommsType,
//            distMap.schedule(),
//            distMap.constructSize(),
//            distMap.subMap(),           // what to send
//            distMap.constructMap(),     // what to receive
//            filmDelta
//        );

        // use solid internal cell as Twall of gas. See if we can make re-radiation stable.
        //kvm, superseded by below  scalarField Twall = nbrIntFld;
        scalarList Twall(patch().size(), 0.0);//kvm
        forAll(*this, i)
        {

            if(0&&filmDelta[i]>1e-6){ //wet
                //fix2
                Twall[i] = max(Tfilm[i],298.15); // film temperature, lagged one time step
                Twall[i] = min(Tfilm[i],378.4); // film temperature, lagged one time step
                //Twall[i] = 444.0; // film temperature
//                if(i%100==0){
//                    /*DEBUG(i);*/
//                    /*DEBUG(Twall[i]);*/
//                    Info << "Twall["<<i<<"] = "<<Twall[i]<<" wet"<<endl;
//                }
            }
            else if(0){ //dry
                Twall[i] = nbrIntFld[i]; // pyrolysis surface temperature (actually, internal temperature)
//                if(i%100==0){
//                    /*DEBUG(i);*/
//                    /*DEBUG(Twall[i]);*/
//                    Info << "Twall["<<i<<"] = "<<Twall[i]<<" dry"<<endl;
//                }
            }
            else{

                scalar filmDeltaWet=0.0002; //m
                scalar filmDeltaDry=0.0000; //m

                scalar ratio=(filmDelta[i]-filmDeltaDry)/(filmDeltaWet-filmDeltaDry);
                ratio=max(ratio,0.0);
                ratio=min(ratio,1.0);

                scalar Twet=min(max(Tfilm[i],298.15),378.4);
                scalar Tdry=nbrIntFld[i];

                Twall[i]=ratio*(Twet-Tdry) + Tdry;
            }
        }
//kvm        const scalarField Twall =
//kvm            (myRadField + myKDelta*intFld + nbrKDelta*nbrIntFld)
//kvm           /(myKDelta + nbrKDelta);

        if (debug)
        {
            scalar Qr = gSum(myRadField*patch().magSf());

            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->dimensionedInternalField().name() << " :" << nl
                << "     radiative heat  [W] : " << Qr << nl
                << "     predicted wallT [K] : " << gAverage(Twall) << nl
                << endl;
        }

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }
    else if (mode == noRad)
    {
        this->refValue() = nbrIntFld;
        this->refGrad() = 0.0;
        this->valueFraction() = nbrKDelta / (nbrKDelta + myKDelta);
    }
    else
    {
        FatalErrorIn
        (
            "turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField"
            "::updateCoeffs()"
        )   << "Illegal mode " << operationModeNames[mode]
            << exit(FatalError);
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qc = gSum(nbrConvFlux*patch().magSf());
        scalar Qr = gSum(nbrRadField*patch().magSf());
        scalar Qt = gSum(nbrTotalFlux*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :" << nl
            << "     convective heat[W] : " << Qc << nl
            << "     radiative heat [W] : " << Qr << nl
            << "     total heat     [W] : " << Qt << nl
            << "     walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField::write
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
    os.writeKeyword("neighbourFieldConvectiveName")<<
        neighbourFieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldConvectiveName")<< fieldConvectiveName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("convectiveScaling") << convectiveScaling_ << token::END_STATEMENT << nl;
//    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadiationCoupledQcWallFunctionMixedSTFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
