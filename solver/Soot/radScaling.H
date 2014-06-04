int size_fPatch = radSc_fPatchName.size();
int size_pPatch = radSc_pPatchName.size();
	
label fLabel[size_fPatch];
scalar fMLR = 0;
for(int i = 0; i < size_fPatch; i++) {
	fLabel[i] = mesh.boundaryMesh().findPatchID(radSc_fPatchName[i]);
    if(fLabel[i]<0){
        FatalErrorIn
               (
                   "radScaling.H"
               )   << "Boundary " << radSc_fPatchName[i] << " not found" << nl
                   << abort(FatalError);
    }
	fMLR = fMLR + gSum(phi.boundaryField()[fLabel[i]]);
}
label pLabel[size_pPatch];
scalar pMLR = 0;
for(int j = 0; j < size_pPatch; j++) {
	pLabel[j] = mesh.boundaryMesh().findPatchID(radSc_pPatchName[j]);
	pMLR = pMLR + gSum(phi.boundaryField()[pLabel[j]]);
}
Info<< "MLR Fuel = " << fMLR << nl;
Info<< "MLR Pyro = " << pMLR << nl;

scalar YsEquiv = 0;
if(fMLR + pMLR != 0) {
	YsEquiv = (radSc_Ys_f*fMLR + radSc_Ys_p*pMLR)/(fMLR + pMLR);
}
Info<< "YsEquiv = " << YsEquiv << nl;

if(YsEquiv > 1e-5) {
	Crs = (radSc_Xr_p + (radSc_Xr_f - radSc_Xr_p)*(YsEquiv - radSc_Ys_p)/(radSc_Ys_f - radSc_Ys_p))/radSc_Xr_f;
}
Info<< "RadScale = " << Crs << nl;
