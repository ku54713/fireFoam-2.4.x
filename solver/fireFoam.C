/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    fireFoam

Description
    Transient solver for fires and turbulent diffusion flames with
    reacting Lagrangian parcels, surface film and pyrolysis modeling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "basicReactingCloud.H"
#include "surfaceFilmModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "solidChemistryModel.H"
#include "pyrolysisModelCollection.H"
#include "MULES.H"
#include "IOmanip.H"

#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "printVersion.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
//    #include "createFieldsSoot.H"
    #include "createPyrolysisModel.H"
    #include "createRadiationModel.H"
    #include "createClouds.H"
    #include "createSurfaceFilmModel.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "readPyrolysisTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << nl;

    while (runTime.run())
    {
        #include "readPISOControls.H"
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"
        // #include "setDeltaT.H"
        #include "readMultivarMULEControls.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << nl;

        parcels.evolve();

        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
            #include "rhoEqn.H"

            // --- Pressure-velocity PIMPLE corrector loop
            for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
            {
                #include "UEqn.H"
                /*if (solveSoot)
                {
                    #include "computeQr.H"
                }*/
                #include "mvConvection.H"
                #include "YhsEqn.H"

                // --- PISO loop
                for (int corr=0; corr<nCorr; corr++)
                {
                    #include "p_rghEqn.H"
                }

                 //turbulence->correct();

                 if (oCorr == nOuterCorr-1)
                 {
                     #include "infoOutput.H"

                     /*if (solveSoot)
                     {
                         #include "computeHp.H"
                     }*/
                 }

                 turbulence->correct();
             }

             rho = thermo.rho();
             kappa = thermo.Cp()*thermo.alpha();
             cp = thermo.Cp();

             if(solvePyrolysisRegion){
                 pyrolysis.evolve();
             }

             runTime.write();
        }
        else
        {
            if(solvePyrolysisRegion){
                pyrolysis.evolve();
            }
            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s";

        if (runTime.outputTime())
        {
            Info<< " +";
        }
        Info<< nl << nl;

    }

    Info<< "End\n" << nl;

    // explicit abort necessary to keep solver from hanging on
    // double-linked list error
    if(Pstream::parRun()){
        Foam::sleep(1);
    }
    abort();

    return(0);
}


// ************************************************************************* //
