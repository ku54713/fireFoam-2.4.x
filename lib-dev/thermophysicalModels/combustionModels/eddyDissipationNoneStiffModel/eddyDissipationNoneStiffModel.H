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
    Foam::combustionModels::eddyDissipationNoneStiffModel

Description
    Eddy dissipation model based on the principle of mixed is burnt. 

SourceFiles
    eddyDissipationNoneStiffModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationNoneStiffModel_H
#define eddyDissipationNoneStiffModel_H
#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationNoneStiffModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationNoneStiffModel
:
    public singleStepCombustion <CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar C_;


    // Private Member Functions

        //- Disallow copy construct
        eddyDissipationNoneStiffModel(const eddyDissipationNoneStiffModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationNoneStiffModel&);


public:

    //- Runtime type information
    TypeName("eddyDissipationNoneStiffModel");


    // Constructors

        //- Construct from components
        eddyDissipationNoneStiffModel
        (
            const word& modelType, const fvMesh& mesh
        );


    //- Destructor
    virtual ~eddyDissipationNoneStiffModel();


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
#   include "eddyDissipationNoneStiffModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
