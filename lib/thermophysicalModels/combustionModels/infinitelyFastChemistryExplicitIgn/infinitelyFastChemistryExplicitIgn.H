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

Class
    Foam::combustionModels::infinitelyFastChemistryExplicitIgn

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    infinitelyFastChemistryExplicitIgn.C

\*---------------------------------------------------------------------------*/

#ifndef infinitelyFastChemistryExplicitIgn_H
#define infinitelyFastChemistryExplicitIgn_H

//#include "combustionModel.H"
//#include "singleStepReactingMixture.H"
//#include "thermoPhysicsTypes.H"
#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class infinitelyFastChemistryExplicitIgn Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class infinitelyFastChemistryExplicitIgn
:
    public singleStepCombustion <CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;

        //- Normalised consumption rate of (fu - fres)
        volScalarField wFuelNorm_;

        //- consumption source term for ignition
        volScalarField wFuelNormIgn_;


    // Private Member Functions

        //- Disallow copy construct
        infinitelyFastChemistryExplicitIgn(const infinitelyFastChemistryExplicitIgn&);

        //- Disallow default bitwise assignment
        void operator=(const infinitelyFastChemistryExplicitIgn&);


public:

    //- Runtime type information
    TypeName("infinitelyFastChemistryExplicitIgn");


    // Constructors

        //- Construct from components
        infinitelyFastChemistryExplicitIgn
        (
            const word& modelType, const fvMesh& mesh
        );


    //- Destructor
    virtual ~infinitelyFastChemistryExplicitIgn();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix, i.e. source term for fuel equation
            virtual tmp<fvScalarMatrix> R(const volScalarField& Y) const;

            //- source term for ignition products
//            virtual tmp<fvScalarMatrix> R2(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //- Return normalised consumption rate of (fu - fres)
            virtual tmp<volScalarField> wFuelNorm() const;

            //- Return ignition consumption source term
            virtual tmp<volScalarField> wFuelNormIgn() const;


    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "infinitelyFastChemistryExplicitIgn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
