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
    Foam::combustionModels::eddyDissipationDiffusionNoneStiffModel

Description
    Eddy dissipation model based on the principle of mixed is burnt. 

SourceFiles
    eddyDissipationDiffusionNoneStiffModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationDiffusionNoneStiffModel_H
#define eddyDissipationDiffusionNoneStiffModel_H
#include "singleStepCombustion.H"

#include "LESModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationDiffusionNoneStiffModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationDiffusionNoneStiffModel
:
    public singleStepCombustion <CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;

        scalar Cd_;

    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- Disallow copy construct
        eddyDissipationDiffusionNoneStiffModel(const eddyDissipationDiffusionNoneStiffModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationDiffusionNoneStiffModel&);

       
public:

    //- Runtime type information
    TypeName("eddyDissipationDiffusionNoneStiffModel");


    // Constructors

        //- Construct from components
        eddyDissipationDiffusionNoneStiffModel
        (
            const word& modelType, const fvMesh& mesh
        );


    //- Destructor
    virtual ~eddyDissipationDiffusionNoneStiffModel();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissipationDiffusionNoneStiffModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
