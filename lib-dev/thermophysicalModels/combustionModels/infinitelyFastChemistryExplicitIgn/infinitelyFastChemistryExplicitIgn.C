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

#include "infinitelyFastChemistryExplicitIgn.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::infinitelyFastChemistryExplicitIgn
(
    const word& modelType, const fvMesh& mesh
)
:
    singleStepCombustion<CombThermoType, ThermoType>(modelType, mesh),
    C_(readScalar(this->coeffs().lookup("C"))),
    wFuelNorm_
    (
        IOobject
        (
            "wFuelNorm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    ),
    wFuelNormIgn_
    (
        IOobject
        (
            "wFuelNormIgn",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::~infinitelyFastChemistryExplicitIgn()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::correct()
{
//    this->singleMixture_.fresCorrect();

    const label fuelI = this->singleMixture_.fuelIndex();

    const volScalarField& YFuel = this->thermo_->composition().Y()[fuelI];

    const dimensionedScalar s = this->singleMixture_.s();

    const dimensionedScalar qFuel_ = this->singleMixture_.qFuel();

    dimensionedScalar nond
    (
            "nond",
            dimMass/(dimLength*pow3(dimTime)),
            1.0
    );

    if (this->thermo_->composition().contains("O2"))
    {
        const volScalarField& YO2 = this->thermo_->composition().Y("O2");
        wFuelNorm_ == this->rho()/(this->mesh_.time().deltaT()*C_)*min(YFuel, YO2/s.value());
        wFuelNormIgn_ == nond/qFuel_;
    }
}


template<class CombThermoType, class ThermoType>
tmp<fvScalarMatrix>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::R(const volScalarField& Y) const
{
    const label specieI = this->thermo_->composition().species()[Y.name()];

    /*const label fNorm = this->singleMixture_.specieProd()[specieI];*/

    const volScalarField fres = this->singleMixture_.fres(specieI);

    const volScalarField wSpecie =
        wFuelNorm_*this->singleMixture_.specieStoichCoeffs()[specieI];
//      / max(fNorm*(Y - fres), SMALL);

    //return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    return wSpecie + fvm::Sp(dimensionedScalar("zero", wSpecie.dimensions(), scalar(0.0)), Y);
}

//template<class CombThermoType, class ThermoType>
//tmp<fvScalarMatrix>
//infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::R2(volScalarField& Y) const
//{
//    const label specieI = this->thermo_->composition().species()[Y.name()];
//
//    const label fNorm = this->singleMixture_.specieProd()[specieI];
//
//    const volScalarField fres = this->singleMixture_.fres(specieI);
//
//    scalar onlyProd(0);
//    if ( Y.name() == "CO2" || Y.name() == "H2O" || Y.name() == "O2" ) {
//        onlyProd = 1.0;
//    }
//    else {
//        onlyProd = 0.0;
//    }
//
//    const volScalarField wSpecie =
//        onlyProd*wFuelNormIgn_*this->singleMixture_.specieStoichCoeffs()[specieI];
////      / max(fNorm*(Y - fres), SMALL);
//
//    //return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
//    return wSpecie + fvm::Sp(dimensionedScalar("zero", wSpecie.dimensions(), scalar(0.0)), Y);
//}
//
template<class CombThermoType, class ThermoType>
tmp<volScalarField>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::dQ() const
{
    const label fuelI = this->singleMixture_.fuelIndex();
    const volScalarField& YFuel = this->thermo_->composition().Y()[fuelI];

    return -this->singleMixture_.qFuel()*(R(YFuel) & YFuel);
}


template<class CombThermoType, class ThermoType>
tmp<volScalarField>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::wFuelNorm() const
{
    return wFuelNorm_;
}

template<class CombThermoType, class ThermoType>
tmp<volScalarField>
infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::wFuelNormIgn() const
{
    return wFuelNormIgn_;
}

template<class CombThermoType, class ThermoType>
bool infinitelyFastChemistryExplicitIgn<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
