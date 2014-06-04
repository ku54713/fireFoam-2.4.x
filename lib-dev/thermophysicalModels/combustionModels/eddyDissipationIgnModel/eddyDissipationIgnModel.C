/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "eddyDissipationIgnModel.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationIgnModel<CombThermoType, ThermoType>::eddyDissipationIgnModel
(
    const word& modelType, const fvMesh& mesh
)
:
    eddyDissipationModel<CombThermoType, ThermoType>(modelType, mesh),
    qIgn_(readScalar(this->coeffs().lookup("qIgn"))),
    chiIgn_(readScalar(this->coeffs().lookup("chiIgn"))),
    tIgnBegin_(readScalar(this->coeffs().lookup("ignStartTime"))),
    tIgnRampUp_(readScalar(this->coeffs().lookup("ignRampUpTime"))),
    tIgnRampDown_(readScalar(this->coeffs().lookup("ignRampDownTime"))),
    tIgnEnd_(readScalar(this->coeffs().lookup("ignEndTime"))),
    wFuelIgn_
     (
         IOobject
         (
             "wFuelIgn",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
     ),
    dQI
     (
         IOobject
         (
             "dQI",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("dQI", dimEnergy/dimTime/pow3(dimLength), 1.1e+06)
     ),
    dQi
     (
         IOobject
         (
             "dQi",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("dQi", dimEnergy/dimTime/pow3(dimLength), 1.1e+06)
     ),
    qrI
     (
         IOobject
         (
             "qrI",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("qrI", dimEnergy/dimTime/pow3(dimLength), 0.0)
     )
{read();}


// * * * * * * * * * * * * * * *lnInclude/combustionModelI.H: * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationIgnModel<CombThermoType, ThermoType>::~eddyDissipationIgnModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void eddyDissipationIgnModel<CombThermoType, ThermoType>::correct()
{
    eddyDissipationModel<CombThermoType, ThermoType>::correct();

    this->wFuelIgn_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        const label cellZoneID = this->mesh().cellZones().findZoneID("igniter");
        const labelList& cellAddr = this->mesh().cellZones()[cellZoneID];

        scalar tRamp_(0); 
        const scalar tCurrent_ = this->db().time().value();

        if(tCurrent_ >= tIgnBegin_ && tCurrent_ <= tIgnRampUp_) {
                tRamp_ = (tCurrent_ - tIgnBegin_)/(tIgnRampUp_ - tIgnBegin_);
        }
        else if(tCurrent_ > tIgnRampUp_ && tCurrent_ < tIgnRampDown_) {
                tRamp_ = 1.0;
        }
        else if(tCurrent_ >= tIgnRampDown_ && tCurrent_ <= tIgnEnd_) {
                tRamp_ = (tIgnEnd_ - tCurrent_)/(tIgnEnd_ - tIgnRampDown_);
        }

        Info << "ign HRR frac = " << tRamp_ << endl;

        volScalarField dQi_ = 0.0*dQI;
        volScalarField qrI_ = 0.0*dQI;

        forAll (cellAddr, cellI)
        {
                dQi_[cellAddr[cellI]] = tRamp_*dQI[cellAddr[cellI]];
                qrI_[cellAddr[cellI]] = chiIgn_*dQi_[cellAddr[cellI]];
        }
        dQi = dQi_;
        qrI = qrI_;

	const dimensionedScalar qFuel_ = this->singleMixture_.qFuel();
        this->wFuelIgn_ == dQi/qFuel_;
    }
}


template<class CombThermoType, class ThermoType>
bool eddyDissipationIgnModel<CombThermoType, ThermoType>::read()
{
	if (this->mesh().cellZones().findZoneID("igniter") == -1)
	{
    	FatalErrorIn
    	(
        	"eddyDissipationIgnModel::New"
    	)   	<< "	Igniter cellZone missing "
        	<< exit(FatalError);
	}

    	const label cellZoneID = this->mesh().cellZones().findZoneID("igniter");
    	const cellZone& zone = this->mesh().cellZones()[cellZoneID];
    	const cellZoneMesh& zoneMesh = zone.zoneMesh();
    	const labelList& cellsZone = zoneMesh[cellZoneID];

    	scalar volCellZone_(0);
    	forAll(cellsZone, cI)
    	{
            	volCellZone_ += this->mesh().V()[cellsZone[cI]];
    	}
    	reduce(volCellZone_,sumOp<scalar>());

    	Info << "Igniter cell volume = " << volCellZone_ << " m^3" << endl;

        dimensionedScalar unitDimScal2("unitDimScal2", dimEnergy/dimTime/pow3(dimLength), 1.0);
        dQI = unitDimScal2*qIgn_/volCellZone_;

        Info << "Igniter HRR/Vol = " << qIgn_/(volCellZone_*1.0e+03) << " kW/m^3" << endl;

	if (dQI[0] > 1.0e+07)
	{
    	FatalErrorIn
    	(
        	"eddyDissipationIgnModel::New"
    	)   	<< "	Igniter HRR/Vol exceeds 10 MW/m^3 "
        	<< exit(FatalError);
	} 

        return true;
}

 template<class CombThermoType, class ThermoType>
 Foam::tmp<Foam::fvScalarMatrix>
 eddyDissipationIgnModel<CombThermoType, ThermoType>::R
 (
     const volScalarField& Y
 ) const
 {
     const label specieI = this->thermo_->composition().species()[Y.name()];

     scalar onlyProd(0);
     if ( Y.name() == "CO2" || Y.name() == "H2O" ) {
         onlyProd = this->singleMixture_.specieStoichCoeffs()[specieI];
     }
     else if ( Y.name() == "O2" ) {
         onlyProd = this->singleMixture_.specieStoichCoeffs()[specieI]-1.0;
     } 
     else {
         onlyProd = 0.0;
     }

     volScalarField wSpecieIgn
     (
         onlyProd*wFuelIgn_ 
     );

     if ( Y.name() == "O2" ) {		// ensuring YO2 doesn't become negative
	wSpecieIgn = max(wSpecieIgn,-this->rho()*Y/this->mesh().time().deltaT());
     }

     return eddyDissipationModel<CombThermoType, ThermoType>::R(Y) + wSpecieIgn + fvm::Sp(0.0*wSpecieIgn, Y);
 }

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationIgnModel< CombThermoType, ThermoType>::Sh() const
{
    return eddyDissipationModel<CombThermoType, ThermoType>::Sh() + dQi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
