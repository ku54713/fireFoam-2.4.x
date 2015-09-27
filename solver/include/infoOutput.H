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

        if (thermo.composition().contains("O2")) 
        {
            ft = (fu*s-O2+YO2Inf)/(s+YO2Inf);
        }
        else
        {
            ft = fu;
        }

        ft.max(0.0);
        ft.min(1.0);

        HRR_fu = ( - mvConvection->interpolate(phi,fu)*phi
                   + (
                         fvc::interpolate
                         (
                             turbulence->alphaEff()*fvc::grad(fu)
                         )
                         & mesh.Sf()
                     )                         
                 ) * qFuel;


        if (mag(g).value() > 0.0) 
        {    
            flameHeight = dimensionedScalar("zero",dimLength,0.0);

            forAll (ft.internalField(),cellI)
            {
                if (ft[cellI] >= ftSt)
                { 
                    flameHeight[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
                }
            }

            outFlameHeight <<  runTime.value() << "  "     
                    << max(flameHeight).value() << "  "
                    <<  nl;
        }


        hc = thermo.hc();  
        h = thermo.he() + hc;

        phiFt = mvConvection->interpolate(phi,ft)*phi
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(ft)*mesh.magSf();

        phiO2 = mvConvection->interpolate(phi,O2)*phi
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(O2)*mesh.magSf();

        phiCO2 = mvConvection->interpolate(phi,CO2)*phi
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(CO2)*mesh.magSf();

        phiHc = phi*mvConvection->interpolate(phi,hc)
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(hc)*mesh.magSf();

        phiHs = phi*mvConvection->interpolate(phi,thermo.he())
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(thermo.he())*mesh.magSf();

        phiH  = phi*mvConvection->interpolate(phi,h)
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(h)*mesh.magSf();


        UT = U * T;
        rhoU = rho * U;
        rhoT = rho * T;
        rhoFu = rho * fu; 

        B = turbulence->R();

        convectiveHeatFlux_L = 
                 - fvc::interpolate(turbulence->alpha())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);

        convectiveHeatFlux_T = 
                 - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);

        //forAll(wallConvectiveHeatFlux.boundaryField(), patchi)
        forAll(QcWallFunction.boundaryField(), patchi)
        {
            //wallConvectiveHeatFlux.boundaryField()[patchi] = convectiveHeatFlux.boundaryField()[patchi];
            
            if (mesh.boundary()[patchi].type() == "mappedWall") 
            {
                //forAll ( wallConvectiveHeatFlux.boundaryField()[patchi],faceI)
                forAll ( QcWallFunction.boundaryField()[patchi],faceI)
                {
                    scalar mlr = - phi.boundaryField()[patchi][faceI]/mesh.boundary()[patchi].magSf()[faceI]*2.5*1000.0; //convert to g/m2/s 
                    if (mlr < 0.1)
                    {
                        //QcWallFunction.boundaryField()[patchi][faceI] = min(max(0,wallConvectiveHeatFlux.boundaryField()[patchi][faceI]),QcThreshold)/QcThreshold*QcFlame;
                        QcWallFunction.boundaryField()[patchi][faceI] = min(max(0,convectiveHeatFlux_L.boundaryField()[patchi][faceI]),QcThreshold)/QcThreshold*QcFlame;
                    }
                    else
                    {
                        //QcWallFunction.boundaryField()[patchi][faceI] = 16000.0 * (mlr/10.0/(Foam::exp(min(500.0,mlr/scalar(10.0)))-scalar(1)));
                        QcWallFunction.boundaryField()[patchi][faceI] = QcFlame * (mlr/10.0/(Foam::exp(mlr/scalar(10.0))-scalar(1)));
                    }
                }
            }
        }

        
/* 
        // alternative way to compute convective heat flux
        forAll(wallConvectiveHeatFlux.boundaryField(), patchi)
        {
            wallConvectiveHeatFlux.boundaryField()[patchi]=
           -turbulence->alphaEff()().boundaryField()[patchi]
           *thermo.Cp()().boundaryField()[patchi]
           *(T.boundaryField()[patchi] - T.boundaryField()[patchi].patchInternalField())
           *mesh.boundary()[patchi].deltaCoeffs();
        }
*/
   
        alphaM = turbulence->alpha();
        muM = turbulence->mu();

        if (!constD) 
        {
            d = turbulence->alpha()/lewisNo/rho;
        }
        else 
        {
            d = DM;
        }
        dSgs = (turbulence->alphaEff()-turbulence->alpha())/rho;
