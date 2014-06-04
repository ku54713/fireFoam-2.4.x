    if( constSR ) {
	C_eta = 0.0;
	srC = 1.0*cSR*C_src;	
    }
    else {
	C_eta = Ceta;
	srC = 0.0*cSR*C_src;	
    }
    strainR = C_eta*pow(turbulence->epsilon()*thermo.rho()/thermo.mu(),0.5) + srC;

    forAll (PsiR, celli) {
	dimensionedScalar SRate2("SRate2",strainR.dimensions(),strainR[celli]);
	dimensionedScalar Hprior2("Hprior2",Hp.dimensions(),Hp[celli]);
	if(SRate2.value() < MinRa0.value()) {
		SRate2.value() = MinRa0.value();
	}
    	if(SRate2.value() >= MinRa0.value() && SRate2.value() <= MaxRa0.value() 
			&& Hprior2.value() >= MinRHp.value() && Hprior2.value() <= MaxRHp.value())
	{ 
		const scalar a0Rsolver = LenRa0.value()*(strainR[celli] - MinRa0.value())/(MaxRa0.value() - MinRa0.value());
    		int i_Ra0 = static_cast<int>(a0Rsolver);
		int i__Ra0 = i_Ra0 + 1;
		const scalar a0Rinterp = (a0Rsolver - scalar(i_Ra0))/(scalar(i__Ra0) - scalar(i_Ra0));
		const scalar HpRsolver = LenRHp.value()*(Hp[celli] - MinRHp.value())/(MaxRHp.value() - MinRHp.value());
    		int j_RHp = static_cast<int>(HpRsolver);
		int j__RHp = j_RHp + 1;
		const scalar HpRinterp = (HpRsolver - scalar(j_RHp))/(scalar(j__RHp) - scalar(j_RHp));
		PsiR[celli] = (1 - a0Rinterp)*(1 - HpRinterp)*psiR[j_RHp][i_Ra0] 
				+ a0Rinterp*(1 - HpRinterp)*psiR[j_RHp][i__Ra0] 
				+ a0Rinterp*HpRinterp*psiR[j__RHp][i__Ra0] 
				+ (1 - a0Rinterp)*HpRinterp*psiR[j__RHp][i_Ra0];
    	}
	else
	{
		PsiR[celli] = 0;
	}
    }

    scalar Crs = 1;
    if( radScaling ) {
	#include "radScaling.H"
    }
    if( Csc > 1.0 && max(thermo.T()).value() > 1000.0) {
    	chiSoot = Csc*Crs*PsiR; 
    }
    else {
    	chiSoot = Csc*Crs*PsiR;	
    }
    chiSoot.min(1);
