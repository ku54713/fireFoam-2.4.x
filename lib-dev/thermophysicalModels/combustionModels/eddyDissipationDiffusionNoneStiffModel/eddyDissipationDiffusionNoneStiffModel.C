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

#include "eddyDissipationDiffusionNoneStiffModel.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::eddyDissipationDiffusionNoneStiffModel
(
    const word& modelType, const fvMesh& mesh
)
:
    singleStepCombustion<CombThermoType, ThermoType>(modelType, mesh),
    C_(readScalar(this->coeffs().lookup("C"))),
    Cd_(readScalar(this->coeffs().lookup("Cd")))
{}


// * * * * * * * * * * * * * * *lnInclude/combustionModelI.H: * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::~eddyDissipationDiffusionNoneStiffModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
//    const compressible::LESModel& sgs = this->mesh().template lookupObject<compressible::LESModel>("LESProperties");
//    return C_*1.048*sqrt(this->turbulence().k())/sgs.delta();
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::rtDiff() const
{
    const compressible::LESModel& sgs = this->mesh().template lookupObject<compressible::LESModel>("LESProperties");

    return Cd_*this->turbulence().alpha()/this->rho()/sqr(sgs.delta());
}


template<class CombThermoType, class ThermoType>
void eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
//        this->singleMixture_.fresCorrect();

        const label fuelI = this->singleMixture_.fuelIndex();

        const volScalarField& YFuel = this->thermo_->composition().Y()[fuelI];

        const dimensionedScalar s = this->singleMixture_.s();

        if (this->thermo_->composition().contains("O2"))
        {
            const volScalarField& YO2 = this->thermo_->composition().Y("O2");

            volScalarField rt = max(rtTurb(),rtDiff()); 

            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                / this->mesh_.time().deltaT()
                * (1 - exp(- this->mesh_.time().deltaT() * rt));
        }
    }
}


template<class CombThermoType, class ThermoType>
bool eddyDissipationDiffusionNoneStiffModel<CombThermoType, ThermoType>::read()
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
