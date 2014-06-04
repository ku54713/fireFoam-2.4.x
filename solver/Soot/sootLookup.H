//#include "IFstream.H"
// Using standard C/C++ syntax to read variables needed for opening PsiSlookup.txt and PsiRlookup.txt files
// Reading these variables from PsiLookup.txt file
  float maxSa0, minSa0, maxSHp, minSHp;
  float maxRa0, minRa0, maxRHp, minRHp;
  int lenSa0, lenSHp;
  int lenRa0, lenRHp;
  if (solveSoot) {
  	FILE * lookupFile;
    fileName psiLookupFile="constant/PsiLookup.txt";
  	lookupFile = fopen (psiLookupFile.c_str(),"r");
  	if(lookupFile==NULL){
        FatalErrorIn
               (
                   "sootLookup"
               )   << "Unable to open file " << "constant/PsiLookup.txt" << nl
                   << abort(FatalError);
  	}
    label ierr(0);
  	ierr = fscanf (lookupFile, "%i %i\n", &lenSa0,&lenSHp);
  	ierr = fscanf (lookupFile, "%f %f %f %f\n", &maxSa0,&minSa0,&maxSHp,&minSHp);
  	ierr = fscanf (lookupFile, "%i %i\n", &lenRa0,&lenRHp);
  	ierr = fscanf (lookupFile, "%f %f %f %f\n", &maxRa0,&minRa0,&maxRHp,&minRHp);
	ierr = fclose (lookupFile);
    if(ierr == EOF ){
            FatalErrorIn("sootLookup.H")
                << "Error reading file" << psiLookupFile
                << abort(FatalError);
    }
  }
  else
  {
	lenSa0 = 0;
	lenSHp = 0;
	maxSa0 = 0;
	minSa0 = 0;
	maxSHp = 0;
	minSHp = 0;
	lenRa0 = 0;
	lenRHp = 0;
	maxRa0 = 0;
	minRa0 = 0;
	maxRHp = 0;
	minRHp = 0;
  }

// Assigning lengths to variables used in reading PsiS and PsiR
  dimensionedScalar LenSa0("LenSa0",dimless,lenSa0);
  dimensionedScalar LenSHp("LenSHp",dimless,lenSHp);
  dimensionedScalar LenRa0("LenRa0",dimless,lenRa0);
  dimensionedScalar LenRHp("LenRHp",dimless,lenRHp);

// Assigning max/min values for strain rate to lookup PsiS values
  dimensionedScalar MaxSa0("MaxSa0",dimless,maxSa0);
  dimensionedScalar MinSa0("MinSa0",dimless,minSa0);

// Assigning max/min values for local heat loss to lookup PsiS values
  dimensionedScalar MaxSHp("MaxSHp",dimless,maxSHp);
  dimensionedScalar MinSHp("MinSHp",dimless,minSHp);

// Assigning max/min values for strain rate to lookup PsiR values
  dimensionedScalar MaxRa0("MaxRa0",dimless,maxRa0);
  dimensionedScalar MinRa0("MinRa0",dimless,minRa0);

// Assigning max/min values for local heat loss to lookup PsiR values
  dimensionedScalar MaxRHp("MaxRHp",dimless,maxRHp);
  dimensionedScalar MinRHp("MinRHp",dimless,minRHp);

  scalar psiS[lenSa0][lenSHp];
  scalar psiR[lenRa0][lenRHp];

  if (solveSoot) {
	// Reading PsiS from PsiSlookup.txt file
	// following OpenFOAM syntax
        fileName fNameS("constant/PsiSlookup.txt");  // This file to be kept in the root of the case directory

        fNameS.expand();

        IFstream PsiSlookupFile(fNameS);

        if (PsiSlookupFile.good())
        {
            for (int i=0; i<lenSa0; i++)
            {
               for (int j=0; j<lenSHp; j++)
               {
                PsiSlookupFile >> psiS[i][j];
               }
            }
        }
        else
        {
            FatalErrorIn("fireMixture")
                << "Cannot read file" << fNameS
                << abort(FatalError);
        }

	// Reading PsiR from PsiRlookup.txt file
	// following OpenFOAM syntax
        fileName fNameR("constant/PsiRlookup.txt"); // This file to be kept in the root of the case directory

        fNameR.expand();

        IFstream PsiRlookupFile(fNameR);

        if (PsiRlookupFile.good())
        {
            for (int k=0; k<lenRa0; k++)
            {
               for (int l=0; l<lenRHp; l++)
               {
                PsiRlookupFile >> psiR[k][l];
		//Info << "psiR = " << psiR[k][l] << " i = " << k << " j = " << l << nl;
               }
            }
        }
        else
        {
            FatalErrorIn("fireMixture")
                << "Cannot read file" << fNameR
                << abort(FatalError);
        }
   }
   else
   {
            for (int i=0; i<lenSa0; i++)
            {
               for (int j=0; j<lenSHp; j++)
               {
                psiS[i][j] = 0;
               }
            }
            for (int i=0; i<lenRa0; i++)
            {
               for (int j=0; j<lenRHp; j++)
               {
                psiR[i][j] = 0;
               }
            }
   }	
