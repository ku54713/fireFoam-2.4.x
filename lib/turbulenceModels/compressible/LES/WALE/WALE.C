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

#include "WALE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WALE, 0);
addToRunTimeSelectionTable(LESModel, WALE, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void WALE::updateSubGridScaleFields()
{
//----------------muSgs is independent of kSgs-------------------------------
//    muSgs_ = ck_*rho()*sqrt(k_)*delta();
//	muSgs_ = rho()*sqr(0.55*delta())*sMut;
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE::WALE
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
//    LESModel(typeName, rho, U, phi, thermoPhysicalModel),
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


    ck_			//ck_ will be useless
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),

    cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cw",
            coeffDict_,
            0.5
        )
    )
{
    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WALE::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

	volSymmTensorField Sij = symm(gradU);
	//volSymmTensorField SS = Sij & Sij;
	volScalarField SuSu = Sij && Sij;
	volTensorField gij = gradU & gradU;
	volTensorField Sd = dev(gij) - skew(gij);
	volScalarField SdSd = Sd && Sd;
	muSgs_ = rho()*sqr(cw_*delta())
		*SdSd*sqrt(SdSd)/(sqr(SuSu)*sqrt(SuSu)
			+ SdSd*sqrt(sqrt(SdSd))
			+ dimensionedScalar("SMALL",dimTime/pow3(sqr(dimTime)),SMALL));

    updateSubGridScaleFields();	//Calculate and update muSgs first;

    volScalarField divU(fvc::div(phi()/fvc::interpolate(rho())));
    volScalarField G(2*muSgs_*(gradU && dev(symm(gradU))));

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho(), k_)
      + fvm::div(phi(), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::SuSp(2.0/3.0*rho()*divU, k_)
      - fvm::Sp(ce_*rho()*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    //bound(k_, k0());
    bound(k_, kMin_);

    //updateSubGridScaleFields();
}


bool WALE::read()
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
