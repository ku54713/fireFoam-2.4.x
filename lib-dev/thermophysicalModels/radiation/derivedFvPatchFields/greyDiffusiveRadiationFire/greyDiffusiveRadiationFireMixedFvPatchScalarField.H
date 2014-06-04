/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::greyDiffusiveRadiationFireMixedFvPatchScalarField

Description
    Radiation temperature specified

SourceFiles
    greyDiffusiveRadiationFireMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef greyDiffusiveRadiationFireMixedFvPatchScalarField_H
#define greyDiffusiveRadiationFireMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
        Class greyDiffusiveRadiationFireMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class greyDiffusiveRadiationFireMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of temperature field
        word TName_;

        //- Emissivity
        scalarField emissivity_;

        //- How to obtain the emissivity
        enum emissivityMode
        {
            LOOKUP,
            SOLIDTHERMO
        };

        static const NamedEnum<emissivityMode, 2> emissivityModeNames_;

        //- Emissivity operation mode
        emissivityMode emissivityMode_;


public:

    //- Runtime type information
    TypeName("greyDiffusiveRadiationFire");


    // Constructors

        //- Construct from patch and internal field
        greyDiffusiveRadiationFireMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        greyDiffusiveRadiationFireMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  greyDiffusiveRadiationFireMixedFvPatchScalarField onto a new patch
        greyDiffusiveRadiationFireMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationFireMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        greyDiffusiveRadiationFireMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationFireMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new greyDiffusiveRadiationFireMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        greyDiffusiveRadiationFireMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationFireMixedFvPatchScalarField&,
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
                new greyDiffusiveRadiationFireMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the temperature field name
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //  adjustment
            word& TName()
            {
                return TName_;
            }

            //- Return the emissivity
            const scalarField& emissivity() const
            {
                return emissivity_;
            }

            //- Return reference to the emissivity to allow adjustment
            scalarField& emissivity()
            {
                return emissivity_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
