dimensionedScalar hcFuel
  (
      "hcFuel",
      dimVelocity*dimVelocity,
      Hu_/Wu_
  );
dimensionedScalar nond
  (
	"nond",
	dimVelocity*dimVelocity,
	1.0
  );

volScalarField hc_ = thermo.hc();
volScalarField hs_ = thermo.he(); 
volScalarField Hp_ = Hp; 

if( constHP ) {
	Hp_ = 0.0;
}
else if( newHP ) {
      Hp_ = (hs_ + hc_ - ft * hcFuel) / (hc_ - ft * hcFuel + 1.0e-4*nond);
}
else {
  	Hp_ = -(hs_ + hc_ - ft * hcFuel) / (nond*(qFuel / (1.0 + stoicRatio)));
  }
Hp = Hp_;

qrflamelet = chiSoot*dQ;
