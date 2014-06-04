/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "WALE2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WALE2, 0);
addToRunTimeSelectionTable(LESModel, WALE2, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void WALE2::updateSubGridScaleFields()
{
//----------------muSgs is independent of kSgs-------------------------------
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE2::WALE2
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName), 
    GenEddyVisc(rho, U, phi, thermoPhysicalModel),

    k_
    (
        IOobject
        ( 
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),


    ck_			//for combustion, keep ck_=sqrt(Cw_)/Ce, then ck_ will be cancel out in combustion model
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.28864  //0.55^2/1.048
        )
    ),

    cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cw",
            coeffDict_,
            0.55
        )
    )
{
    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WALE2::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

	volSymmTensorField	Sij	=	symm(gradU);
	volScalarField	SuSu	=	Sij && Sij;
	volTensorField	gij	=	gradU & gradU;
	volTensorField	Sd	=	dev(gij) - skew(gij);
	volScalarField	SdSd	=	Sd && Sd;
			muSgs_	=	rho()*sqr(cw_*delta())
					*SdSd*sqrt(SdSd)/(sqr(SuSu)*sqrt(SuSu)
					+ SdSd*sqrt(sqrt(SdSd))
					+ dimensionedScalar("SMALL",dimTime/pow3(sqr(dimTime)),SMALL));

    updateSubGridScaleFields();	//Calculate and update muSgs first;

    k_=sqr(muSgs_/(ck_*rho()*delta()));
    k_.correctBoundaryConditions();
}


bool WALE2::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        cw_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
