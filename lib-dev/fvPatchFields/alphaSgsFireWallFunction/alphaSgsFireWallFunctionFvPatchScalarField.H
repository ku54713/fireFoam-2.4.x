/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Class
    alphaSgsFireWallFunctionFvPatchScalarField

Description
    Buoyant Thermal wall function for turbulent thermal diffusivity based on
    "Asymptotic analysis of the near-wall region of turbulence natural
    convection flows"
    Holling and Herwig. J. Fluid Mech (2005), vol 541, pp 383-397

SourceFiles
    alphaSgsJayatillekeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaSgsFireWallFunctionFvPatchScalarField_H
#define alphaSgsFireWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
     Class alphaSgsFireWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphaSgsFireWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

/*
        static scalar maxExp_;
        static scalar tolerance_;
        static label maxIters_;
*/

        //- Thermal expansion coefficient
        scalar QcFlame_;

        //- gravity magnitude
        scalar QcThreshold_;

    // Private member functions

        //- Check the type of the patch
//        void checkType();

        //- Read
        //void read();
public:

    //- Runtime type information
    TypeName("alphaSgsFireWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphaSgsFireWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaSgsFireWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //  alphaSgsFireWallFunctionFvPatchScalarField
        //  onto a new patch
        alphaSgsFireWallFunctionFvPatchScalarField
        (
            const alphaSgsFireWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaSgsFireWallFunctionFvPatchScalarField
        (
            const alphaSgsFireWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaSgsFireWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphaSgsFireWallFunctionFvPatchScalarField
        (
            const alphaSgsFireWallFunctionFvPatchScalarField&,
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
                new alphaSgsFireWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::blocking
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
