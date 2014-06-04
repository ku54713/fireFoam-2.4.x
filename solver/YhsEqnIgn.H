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
    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
            fvScalarMatrix R = combustion->R(Yi);
            dimensionedScalar nond
            (
                "nond",
                (dimLength*pow3(dimTime))/dimMass,
                1.0
            );
//            fvScalarMatrix R2 = nond*dQi*combustion->R2(Yi);

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
                  - fvm::laplacian(turbulence->alphaEff(), Yi)
                 == 
                    parcels.Srho(i)
                    + surfaceFilm.Srho(i)
//                    + R2
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
                  - fvm::laplacian(turbulence->alphaEff(), Yi)
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
                    parcels.Srho(i)
                    + surfaceFilm.Srho(i)
//                    + R2
                    + R,
                    mesh.solver("Yi")
                );

            }

            //Info << "max of " << Y[i].name() << " = " << max(Y[i]).value() << nl;
            //Info << "min of " << Y[i].name() << " = " << min(Y[i]).value() << nl;

            //Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
        Info << "max of " << Y[i].name() << " = " << max(Y[i]).value() << nl;
        Info << "min of " << Y[i].name() << " = " << min(Y[i]).value() << nl;
    }

    Y[inertIndex] = scalar(1) - Yt;
    //Y[inertIndex].max(0.0);

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
          + dQ                    
          + radiation->Shs(thermo)
      + dQi
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
          + dQ                    
          + radiation->Shs(thermo)
          + dQi
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
          + dQ                    
          + radiation->Shs(thermo)
          + dQi
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
