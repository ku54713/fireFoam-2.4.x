tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
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
            fvScalarMatrix R(combustion->R(Yi));

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              // FIXME: species uses dEff and energy uses alphaEff
              - fvm::laplacian(dEff, Yi)
              ==
                parcels.SYi(i, Yi)
              + surfaceFilm.Srho(i)
              + combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

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
             << setw(12) << setprecision(5) << gMin(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gAverage(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gMax(Y[i]) << nl;
    }

    radiation->correct();

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        combustion->Sh()
      + radiation->Sh(thermo)
      + parcels.Sh(he)
      + surfaceFilm.Sh()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
