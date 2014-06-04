/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::reactingOneDim21CharOxi

Description
    Reacting, 1-D pyrolysis model

SourceFiles
    reactingOneDim21CharOxi.C

\*---------------------------------------------------------------------------*/

#ifndef reactingOneDim21CharOxi_H
#define reactingOneDim21CharOxi_H

#include "reactingOneDim21.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class reactingOneDim21CharOxi Declaration
\*---------------------------------------------------------------------------*/

class reactingOneDim21CharOxi
:
    public reactingOneDim21
{
private:


protected:

    // Protected data

            //- Total char volume fraction [-], only computed at patch
            volScalarField Xchar_;

            //- Total char volume fraction [kg], only computed at patch
            volScalarField mChar_;

            //- Total char volume fraction [kg], only computed at patch
            volScalarField mCharBurnt_;

            //- Heat release due to char oxidation [J/s/m3]
            volScalarField charOxiSh_;

            //- Total O2 mass flux from the primary region [kg/s]
            volScalarField phiO2_;

            //- Total CO2 mass flux to the primary region  [kg/s]
            volScalarField phiCO2_;
 
    // Protected member functions

        //- Update submodels
        virtual void updateFields();

        //- Update char oxidation rate
        virtual void updateCharOxi();

        // Equations

            //- Solve energy
            virtual void solveEnergy();

public:

    //- Runtime type information
    TypeName("reactingOneDim21CharOxi");


    // Constructors

        //- Construct from type name and mesh
        reactingOneDim21CharOxi
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name, mesh and dictionary
        reactingOneDim21CharOxi
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );


    //- Destructor
    virtual ~reactingOneDim21CharOxi();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "reactingOneDim21CharOxiI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
