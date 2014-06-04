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
    Foam::combustionModels::eddyDissipationIgnModel

Description
    Volumetric igniter model based on the eddy dissipation model.

SourceFiles
    eddyDissipationIgnModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationIgnModel_H
#define eddyDissipationIgnModel_H
#include "eddyDissipationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationIgnModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationIgnModel
:
    public eddyDissipationModel <CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
        scalar qIgn_;
        scalar chiIgn_;
    	scalar tIgnBegin_;
    	scalar tIgnRampUp_;
    	scalar tIgnRampDown_;
    	scalar tIgnEnd_;


    // Private Member Functions

        //- Disallow copy construct
        eddyDissipationIgnModel(const eddyDissipationIgnModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationIgnModel&);


protected:

	// Protected data

		//- Igniter fuel consumption rate
		volScalarField wFuelIgn_;

		//- Temporary dQ field
		volScalarField dQI;

public:

    //- Runtime type information
    TypeName("eddyDissipationIgnModel");


    // Constructors

        //- Construct from components
        eddyDissipationIgnModel
        (
            const word& modelType, const fvMesh& mesh
        );


    //- Destructor
    virtual ~eddyDissipationIgnModel();


    // Member Functions

        // Evolution

	    //- Fuel consumption rate matrix
	    virtual tmp<fvScalarMatrix> R(const volScalarField& Y) const;

            //- Correct combustion rate
            virtual void correct();

	    //- Igniter heat release rate
	    volScalarField dQi;

	    //- Igniter radiant power
	    volScalarField qrI;

	    //- Sensible enthalpy source term
	    virtual tmp<volScalarField> Sh() const;

    // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissipationIgnModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
