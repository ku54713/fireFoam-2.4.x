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

#include "SIGMA2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SIGMA2, 0);
addToRunTimeSelectionTable(LESModel, SIGMA2, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void SIGMA2::updateSubGridScaleFields()
{
//----------------muSgs is independent of kSgs-------------------------------

    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SIGMA2::SIGMA2
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


    ck_			//for combustion, keep ck_=sqrt(Cx_)/cw, then ck_ will be cancel out in combustion model
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            1.739   //1.35^2/1.048
        )
    ),

    cx_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cx",
            coeffDict_,
            1.35
        )
    )
{
    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SIGMA2::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);
    
    volSymmTensorField	Gij	=	symm(gradU.T() & gradU);
    volScalarField	L1	=	tr(Gij);
    volScalarField	L2	=	0.5*(sqr(L1) - tr(Gij & Gij));
    dimensionedScalar			sMin4("sMin4",L2.dimensions(),1);
    volScalarField	L3	=	sMin4*det(Gij);

    volScalarField	A1	=	sqr(L1)/9.0 - L2/3.0 + dimensionedScalar("SMALL",L2.dimensions(),SMALL);
    volScalarField	A2	=	L1*sqr(L1)/27.0 - L1*L2/6.0 + L3/2.0;
    volScalarField	R12	=	A2/( pow(A1,1.5) + dimensionedScalar("SMALL",A2.dimensions(),SMALL) );
    volScalarField	A3	=	acos(0.9999*R12)/3.0;

    volScalarField	ss1	=	L1/3.0 + 2.0*sqrt(A1)*cos(A3);
    volScalarField	ss2	=	L1/3.0 - 2.0*sqrt(A1)*cos(3.1415926/3.0+A3);
    volScalarField	ss3	=	L1/3.0 - 2.0*sqrt(A1)*cos(3.1415926/3.0-A3);

    volScalarField	S1	=	sqrt(mag(ss1));
    volScalarField	S2	=	sqrt(mag(ss2));
    volScalarField	S3	=	sqrt(mag(ss3));    

    volScalarField muSgs_real	=	rho()*sqr(cx_*delta())
    					*(S3*(S1-S2)*(S2-S3))
    					/(sqr(S1) + dimensionedScalar("SMALL",sqr(S1.dimensions()),SMALL));
    		
    if(min(muSgs_real).value()<0)
    {
    	Info<<"Problems about muSgs, negtive value:----> "<<min(muSgs_real).value()<<endl;
    }

    muSgs_=mag(muSgs_real);
    updateSubGridScaleFields();	//Calculate and update muSgs first;

    k_=sqr(muSgs_/(ck_*rho()*delta()));
    k_.correctBoundaryConditions();
}


bool SIGMA2::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        cx_.readIfPresent(coeffDict());

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
