/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "thermoSingleLayer.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "massAbsorptionModel.H"
#include "filmRadiationModel.H"
#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchI)
    {
        if (bTypes[patchI] == mappedFieldFvPatchField<scalar>::typeName)
        {
            bTypes[patchI] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // no additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::resetPrimaryRegionSourceTerms()" << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar("zero", hsSp_.dimensions(), 0.0);
}


void thermoSingleLayer::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    mu_ == filmThermo_->mu(); // OpenFOAM version mission this, bug, kvm
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();

    forAll(rho_, cellI)
    {
        const scalar Ts = max(min(Ts_[cellI],Tmax_),Tmin_);
        const scalar p = pPrimary_[cellI];
    }

    rho_.correctBoundaryConditions();
    mu_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
    Cp_.correctBoundaryConditions();
    kappa_.correctBoundaryConditions();

}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    forAll(T_.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchI];
        if (isA<mappedFieldFvPatchField<scalar> >(Tp))
        {
            hs_.boundaryField()[patchI] == hs(Tp, patchI);
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
        if(pyrCoupled_){
            // UIndirectList<scalar>(Tw_, pp.faceCells()) =
            // 		//T_.boundaryField()[patchI];
            // 		pyrTemperaturePtr_->boundaryField()[patchI];
            // TODO: how do I do this with the new coupling?

            scalarList Tpyr(pp.faceCells().size(), 0.0);


            typedef regionModels::pyrolysisModels::pyrolysisModel
                pyrolysisModelType;

            const regionModels::regionModel& pyrolysisRegion =
            db().time().lookupObject<regionModels::regionModel>
                (
                    "pyrolysisProperties"
                    );

            const pyrolysisModelType& pyrolysisModel =
                dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
            
            pyrolysisModelType& pyrolysis =
                const_cast<pyrolysisModelType&>(pyrolysisModel);
            
            // internal cell tempertaure must be used for stability
            Tpyr = 
                mapRegionPatchInternalField<scalar>
                (
                    pyrolysis,
                    "T",
                    patchI,
                    true
                    );
            
            UIndirectList<scalar>(Tw_,pp.faceCells()) = 
                Tpyr;
        }
        else if(1){
            // compute Tw based on 0D lumped capacitance model
            forAll(Tw_,i){
                Tw_[i]=qFilmToWall_[i]*time_.deltaTValue()/2702.0/949.0/0.0012 + Tw_[i];
            }
            Info << "max Tw " << gMax(Tw_) << endl;
        }
        else{
            UIndirectList<scalar>(Tw_, pp.faceCells()) =
            		T_.boundaryField()[patchI];
        }
    }
    Tw_.correctBoundaryConditions();

    // Update heat transfer to wall (used in film/pyrolysis coupling)
    // heat flow out of film is positive
    qFilmToWall_ = htcw_->h()*(T_ - Tw_);

    // Info << "Tp " << T_[10] << endl;
    // Info << "Twp " << Tw_[10] << endl;
    // Info << "htc " << htcw_->h() << endl;

    forAll(qFilmToWall_,i){
        if(delta_[i]<1e-8){
            qFilmToWall_[i]=0.0;
        }
    }

    //scalarField qw(htcw_->h()*(T_-Tw_));
    //scalarField qs(htcs_->h()*(T_-TPrimary_));
    //DEBUG(qw[2]);
    //DEBUG(qs[2]);

    qFilmToWall_.correctBoundaryConditions();

    // Update heat transfer from gas phase (used in diagnostics)
    // heat flow out of film is positive
    qGasToFilm_ = htcs_->h()*(T_ - TPrimary_);
    forAll(qGasToFilm_,i){
        if(delta_[i]<1e-8){
            qGasToFilm_[i]=0.0;
        }
    }
    qGasToFilm_.correctBoundaryConditions();

    // Update film surface temperature
    // TODO: Ts estimation should go here

    const scalarField d2k(delta_/2.0/kappa_);
    const scalarField Shs(radiation_->ShsConst());
    const scalarField hInf(htcs_->h());

    Ts_.internalField() = 
        (T_+d2k*(hInf*TPrimary_+Shs))
        /(1+d2k*hInf);

    // Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::transferPrimaryRegionThermoFields()" << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::transferPrimaryRegionSourceFields()" << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    // Convert accummulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimary_.boundaryField(), patchI)
    {
        const scalarField& priMagSf =
            primaryMesh().magSf().boundaryField()[patchI];

        hsSpPrimary_.boundaryField()[patchI] /= priMagSf*deltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();

    // Apply enthalpy source as difference between incoming and actual states
    hsSp_ -= rhoSp_*hs_;
}


void thermoSingleLayer::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}


void thermoSingleLayer::updateSubmodels()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::updateSubmodels()" << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update vaporization
    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassPCTrans_,
        primaryEnergyPCTrans_
    );

    //for diagnostics
    primaryMassPCTrans_.correctBoundaryConditions();
    primaryEnergyPCTrans_.correctBoundaryConditions();

    // Update massAbsorption
    massAbsorption_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        massAbs_,
        energyAbs_
    );

    //for diagnostics
    massAbs_.correctBoundaryConditions();
    energyAbs_.correctBoundaryConditions();

    // Update radiation
    radiation_->correct();

    // Update kinematic sub-models
    kinematicSingleLayer::updateSubmodels();

    // Update source fields (phase change)
    hsSp_ += primaryEnergyPCTrans_/magSf()/time().deltaT();
    rhoSp_ += primaryMassPCTrans_/magSf()/time().deltaT();

    // Update source fields (mass absorption)
    rhoSp_ += massAbs_/magSf()/time().deltaT();

    // Vapour recoil pressure (can become unstable for wild oscillations in vaporization rate)
    // pSp_ -= sqr(primaryMassPCTrans_/magSf()/time_.deltaT())/2.0/rhoPrimary_;
    // Info << "vaporRecoilPressure " << gMin(pSp_) << " " << gAverage(pSp_) << " " << gMax(pSp_) << nl;
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{
    dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);

   return
       (
        - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
        - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
        );
}


void thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::solveEnergy()" << endl;
    }

    updateSurfaceTemperatures();

    solve
        (
         fvm::ddt(deltaRho_, hs_)
         + fvm::div(phi_, hs_)
         ==
         // is vaporization energy accounted for twice?  hsSp_ and rhoSp_*hs_?
         - hsSp_
         + q(hs_)
         + radiation_->Shs()
         //      - fvm::SuSp(rhoSp_, hs_)
         - rhoSp_*hs_
         );
    
    correctThermoFields();

    // evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    Cp_
    (
        IOobject
        (
            "Cp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "kappa",
            dimEnergy/dimTime/dimLength/dimTemperature,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    hs_
    (
        IOobject
        (
            "hsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0),
        hsBoundaryTypes()
    ),
    qGasToFilm_
    (
     IOobject
     (
      "qGasToFilm",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qFilmToWall_
    (
     IOobject
     (
      "qFilmToWall",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),

    massAbs_
    (
        IOobject
        (
            "massAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    energyAbs_
    (
        IOobject
        (
            "energyMassAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    primaryMassPCTrans_
    (
        IOobject
        (
            "primaryMassPCTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    primaryEnergyPCTrans_
    (
        IOobject
        (
            "primaryEnergyPCTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaWet_(readScalar(coeffs_.lookup("deltaWet"))),
    hydrophilic_(readBool(coeffs_.lookup("hydrophilic"))),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // must have same name as hSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", hsSp_.dimensions(), 0.0)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    viscosity_(filmViscosityModel::New(*this, coeffs(), mu_)),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    massAbsorption_(massAbsorptionModel::New(*this, coeffs())),
    radiation_(filmRadiationModel::New(*this, coeffs())),
    Tmin_(-VGREAT),
    Tmax_(VGREAT)
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar("zero", dimless, 0.0),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.lookup("hydrophilicDryScale") >> hydrophilicDryScale_;
        coeffs_.lookup("hydrophilicWetScale") >> hydrophilicWetScale_;
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();
        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);
        // delta_.correctBoundaryConditions();
        // U_.correctBoundaryConditions();
        deltaRho_ == delta_*rho_;
        phi_ = fvc::interpolate(deltaRho_*U_) & regionMesh().Sf();
        // Info << phi_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchI,
    const label faceI,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchI,
        faceI,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryField()[patchI][faceI] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::preEvolveRegion()" << endl;
    }

//    correctHsForMappedT();

    kinematicSingleLayer::preEvolveRegion();

    // Update phase change
    primaryMassPCTrans_ == dimensionedScalar("zero", dimMass, 0.0);
    primaryEnergyPCTrans_ == dimensionedScalar("zero", dimEnergy, 0.0);
    // Update mass absorption
    massAbs_ == dimensionedScalar("zero", dimMass, 0.0);
    energyAbs_ == dimensionedScalar("zero", dimEnergy, 0.0);
}


void thermoSingleLayer::evolveRegion()
{
    if (debug)
    {
        Info<< "thermoSingleLayer::evolveRegion()" << endl;
    }

    correctAlpha();

    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update temperature using latest hs_
    T_ == T(hs_);

#include "diagnostics.H"

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Update film wall and surface temperatures
    // updateSurfaceTemperatures();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


tmp<volScalarField> thermoSingleLayer::massAbs() const
{
    return massAbs_;
}


tmp<volScalarField> thermoSingleLayer::primaryMassTrans() const
{
    return primaryMassPCTrans_;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = T_.internalField();

    Info<< indent << "min/mean/max(T)    = " << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", " << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
    massAbsorption_->info(Info);
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& Srho = tSrho();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchI = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassPCTrans_.boundaryField()[filmPatchI];

        toPrimary(filmPatchI, patchMass);

        const label primaryPatchI = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchI = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassPCTrans_.boundaryField()[filmPatchI];

            toPrimary(filmPatchI, patchMass);

            const label primaryPatchI = primaryPatchIDs()[i];
            const unallocLabelList& cells =
                primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Sh() const
{
    tmp<DimensionedField<scalar, volMesh> > tSh
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
/*
    phase change energy fed back into the film...

    scalarField& Sh = tSh();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchI = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyPCTrans_.boundaryField()[filmPatchI];

        toPrimary(filmPatchI, patchEnergy);

        const label primaryPatchI = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }
*/
    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
