/*tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_hs)")
    )
);*/
{
    volScalarField dEff
    (
        "dEff",
        turbulence->alphaEff() 
    );
   
    if (!constD) 
    {
        dEff -= turbulence->alpha()*(1-1.0/lewisNo);
    }
    else
    {
        dEff = dEff - turbulence->alpha() + DM*rho;
    }

    combustion->correct();

    // dQ is only for diagnostic purposes.
    // TODO: set dQ = combustoin->dQ(); once OpenCFD updates 
    //       greyAbsorptionEmission model in OpenFOAM source
    dQ = combustion->Sh();

    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
            fvScalarMatrix R = combustion->R(Yi);

            if (isMULES)
            {
                MULES::explicitSolve
                (
                        rho,
                        Yi,
                        phi,
                        fvc::flux
                        (
                        		phi,
                        		Yi
                        )(),
                        zeroField(),
                        zeroField(),
                        1, 0
                );

                solve
                (
                        fvm::ddt(rho, Yi) - fvc::ddt(rho, Yi)
                      //- fvm::laplacian(turbulence->alphaEff(), Yi)
                        - fvm::laplacian(dEff, Yi)
                        == 
                        //parcels.Srho(i)
                        parcels.SYi(i, Yi)
                        + surfaceFilm.Srho(i)
                        + R,
                        mesh.solver("Yi")
                );
            }
            else
            {
                fvScalarMatrix YiEqn
                (
                        fvm::ddt(rho, Yi)
                        + mvConvection->fvmDiv(phi, Yi)
                      //- fvm::laplacian(turbulence->alphaEff(), Yi)
                        - fvm::laplacian(dEff, Yi)
                );

                if (oCorr == nOuterCorr-1)
                {
                    YiEqn.relax(1);
                }
                else
                {
                    YiEqn.relax();
                }

                solve
                (
                        YiEqn
                        ==
                        //parcels.Srho(i)
                        parcels.SYi(i, Yi)
                        + surfaceFilm.Srho(i)
                        + R,
                        mesh.solver("Yi")
                );

            }

            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] == scalar(1) - Yt;
    //Y[inertIndex].max(0.0);

    forAll(Y, i){
        Info << setw(8) << Y[i].name() << "\tmin/ave/max "<< " =\t" 
             << setw(12) << gMin(Y[i]) << "\t" 
             << setw(12) << gAverage(Y[i]) << "\t" 
             << setw(12) << gMax(Y[i]) << nl;
    }

    radiation->correct();

    if (!isMULES)
    {
        fvScalarMatrix hsEqn
        (
                fvm::ddt(rho, hs)
                + mvConvection->fvmDiv(phi, hs)
                - fvm::laplacian(turbulence->alphaEff(), hs)
                ==
                DpDt
                + parcels.Sh(hs)
                + surfaceFilm.Sh()
                // TODO: revert to using combustion->Sh() once 
                //       greyMeanAbsorptionEmission updated
                // + combustion->Sh()
                + XEff*dQ
                //+ radiation->Shs(thermo)
                + radiation->Sh(thermo)
        );

        if (oCorr == nOuterCorr-1)
        {
            hsEqn.relax(1);
        }
        else
        {
            hsEqn.relax();
        }

        hsEqn.solve();
    }
    else if (isHsMULES)
    {
        MULES::explicitSolve
        (
                rho,
                hs,
                phi,
                fvc::flux
                (
                        phi,
                        hs
                )(),
                zeroField(),
                zeroField(),
                //1, 0
                1e10, 0
        );

        solve
        (
                fvm::ddt(rho, hs) - fvc::ddt(rho, hs)
                - fvm::laplacian(turbulence->alphaEff(), hs)
                ==
                DpDt
                + parcels.Sh(hs)
                + surfaceFilm.Sh()
                // + combustion->Sh()
                + XEff*dQ
                //+ radiation->Shs(thermo)
                + radiation->Sh(thermo)
        );
    }
    else
    {
        fvScalarMatrix hsEqn
        (
                fvm::ddt(rho, hs)
                + fvm::div(phi, hs)
                - fvm::laplacian(turbulence->alphaEff(), hs)
                ==
                DpDt
                + parcels.Sh(hs)
                + surfaceFilm.Sh()
                // + combustion->Sh()
                + XEff*dQ
                //+ radiation->Shs(thermo)
                + radiation->Sh(thermo)
        );

        if (oCorr == nOuterCorr-1)
        {
            hsEqn.relax(1);
        }
        else
        {
            hsEqn.relax();
        }

        hsEqn.solve();
    }

    thermo.correct();
}
