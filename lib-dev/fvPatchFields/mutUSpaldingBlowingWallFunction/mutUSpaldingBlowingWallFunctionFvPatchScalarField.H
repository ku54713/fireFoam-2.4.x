/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::compressible::mutUSpaldingBlowingWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    This boundary condition provides a turbulent viscosity condition when
    using wall functions for rough walls, based on velocity,  using  Spalding's
    law to give a continuous nut profile to the wall (y+ = 0)

        \f[
            y^+ = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
                - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
        \f]

    where
    \vartable
        y^+     | non-dimensional position
        u^+     | non-dimensional velocity
        \kappa  | Von Karman constant
    \endvartable


    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            mutUSpaldingBlowingWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::compressible::mutWallFunctionFvPatchScalarField

SourceFiles
    mutUSpaldingBlowingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutUSpaldingBlowingWallFunctionFvPatchScalarField_H
#define compressibleMutUSpaldingBlowingWallFunctionFvPatchScalarField_H

//#include "mutWallFunctionFvPatchScalarField.H"
#include "mutUSpaldingWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class mutUSpaldingBlowingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutUSpaldingBlowingWallFunctionFvPatchScalarField
:
    public mutUSpaldingWallFunctionFvPatchScalarField
{

     // Private data

        //- h over Cp (default = 10), unit (W/m^2/K)/(J/g/K)=g/s/m^2
        scalar hOCp_;

        //- ratio between phi_pyrolysate and phi_fuel (phiGas/phi), default = 2.5 for corrugated over C3H8
        scalar rPhi_;

protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;
/*
        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;
*/

public:

    //- Runtime type information
    TypeName("mutUSpaldingBlowingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutUSpaldingBlowingWallFunctionFvPatchScalarField
        //  onto a new patch
        mutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const mutUSpaldingBlowingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const mutUSpaldingBlowingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutUSpaldingBlowingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const mutUSpaldingBlowingWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new mutUSpaldingBlowingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
/*
        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;

*/

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************ //
