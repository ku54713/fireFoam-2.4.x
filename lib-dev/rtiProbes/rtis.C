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

\*---------------------------------------------------------------------------*/

#include "rtis.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include <dirent.h>
#include "IFstream.H"
#include "IStringStream.H"
#include "interpolationTable.H"

#include <istream>
#include <sys/stat.h>
#include <vector>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(rtis, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rtis::findElements(const fvMesh& mesh)
{
    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    rtiValues_.clear();
    rtiValues_.setSize(size());

    dTeOld_.clear();
    dTeOld_.setSize(size(),0.0);

    forAll(*this, rtiI)
    {
        const vector& location = operator[](rtiI);

        const label cellI = mesh.findCell(location);

        elementList_[rtiI] = cellI;

        if (cellI != -1)
        {
            const labelList& cellFaces = mesh.cells()[cellI];
            const vector& cellCentre = mesh.cellCentres()[cellI];
            scalar minDistance = GREAT;
            label minFaceID = -1;
            forAll (cellFaces, i)
            {
                label faceI = cellFaces[i];
                vector dist = mesh.faceCentres()[faceI] - cellCentre;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = faceI;
                }
            }
            faceList_[rtiI] = minFaceID;
        }
        else
        {
            faceList_[rtiI] = -1;
        }

        if (debug && (elementList_[rtiI] != -1 || faceList_[rtiI] != -1))
        {
            Pout<< "rtis : found point " << location
                << " in cell " << elementList_[rtiI]
                << " and face " << faceList_[rtiI] << endl;
        }
    }


    // Check if all rtis have been found.
    forAll(elementList_, rtiI)
    {
        const vector& location = operator[](rtiI);
        label cellI = elementList_[rtiI];
        label faceI =  faceList_[rtiI];

        // Check at least one processor with cell.
        reduce(cellI, maxOp<label>());
        reduce(faceI, maxOp<label>());

        if (cellI == -1)
        {
            if (Pstream::master())
            {
                WarningIn("findElements::findElements(const fvMesh&)")
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (faceI == -1)
        {
            if (Pstream::master())
            {
                WarningIn("rtis::findElements(const fvMesh&)")
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[rtiI] != -1 && elementList_[rtiI] != cellI)
            {
                WarningIn("rtis::findElements(const fvMesh&)")
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[rtiI]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << cellI << " on some other domain."
                    << endl
                    << "This might happen if the rti location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[rtiI] != -1 && faceList_[rtiI] != faceI)
            {
                WarningIn("rtis::findElements(const fvMesh&)")
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[rtiI]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << faceI << " on some other domain."
                    << endl
                    << "This might happen if the rti location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::rtis::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        if (debug)
        {
            Info<< "Probing fields:" << currentFields << nl
                << "Probing locations:" << *this << nl
                << endl;
        }


        fileName rtiDir;
        fileName rtiSubDir = name_;

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            rtiSubDir = rtiSubDir/mesh_.name();
        }
        rtiSubDir = rtiSubDir/mesh_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            rtiDir = mesh_.time().path()/".."/rtiSubDir;
        }
        else
        {
            rtiDir = mesh_.time().path()/rtiSubDir;
        }

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, rtiFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                if (debug)
                {
                    Info<< "close rti stream: " << iter()->name() << endl;
                }

                delete rtiFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(rtiDir);

            OFstream* sPtr = new OFstream(rtiDir/fieldName);

            if (debug)
            {
                Info<< "open rti stream: " << sPtr->name() << endl;
            }

            rtiFilePtrs_.insert(fieldName, sPtr);

            unsigned int w = IOstream::defaultPrecision() + 7;

            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                     << vector::componentNames[cmpt];

                forAll(*this, rtiI)
                {
                    *sPtr<< ' ' << setw(w) << operator[](rtiI)[cmpt];
                }
                *sPtr << endl;
            }

            *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                 << "Time" << endl;
        }
    }

    prepareRTI();

    return nFields;
}


void Foam::rtis::prepareRTI()
{
    // adjust file streams
    if (Pstream::master())
    {

        fileName rtiDir;
        fileName rtiSubDir = name_;

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            rtiSubDir = rtiSubDir/mesh_.name();
        }
        rtiSubDir = rtiSubDir/mesh_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            rtiDir = mesh_.time().path()/".."/rtiSubDir;
        }
        else
        {
            rtiDir = mesh_.time().path()/rtiSubDir;
        }

        const word& fieldName = "RTI";

        // Create directory if does not exist.
        mkDir(rtiDir);
        
        OFstream* sPtr = new OFstream(rtiDir/fieldName);

        if (debug)
        {
            Info<< "open rti stream: " << sPtr->name() << endl;
        }
	
        rtiFilePtrs_.insert(fieldName, sPtr);

        unsigned int w = IOstream::defaultPrecision() + 7;

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                 << vector::componentNames[cmpt];
	    
            forAll(*this, rtiI)
            {
                *sPtr<< ' ' << setw(w) << operator[](rtiI)[cmpt];
            }
            *sPtr << endl;
        }
	
        *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
             << "Time" << endl;
    }
    
}

void Foam::rtis::writeRTI()
{
    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *rtiFilePtrs_["RTI"];

        os  << setw(w) << mesh_.time().value();

        forAll(rtiValues_, rtiI)
        {
            os  << ' ' << setw(w) << rtiValues_[rtiI];
        }
        os  << endl;
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rtis::rtis
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    pointField(0),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    cleanRestart_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rtis::~rtis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rtis::execute()
{

    
    const volScalarField& T = mesh_.lookupObject<volScalarField> ("T");
    const volVectorField& U = mesh_.lookupObject<volVectorField> ("U");
    forAll(*this, rtiI)
    {
        label cellI = elementList_[rtiI];
        if(cellI!=-1)
        {
            rtiValues_[rtiI]=calculateLinkTemperature(rtiI,T[cellI],mag(U[cellI]));
        }
        else
        {
            rtiValues_[rtiI]=-1.0;
        }
    }
    // bring local values back to each processor
    forAll(*this,rtiI)
    {
        reduce(rtiValues_[rtiI], maxOp<scalar>());
    }
    // writeRTI();
    // Do nothing - only valid on write
}

Foam::scalar Foam::rtis::calculateLinkTemperature(label rtiI, scalar Tgas, scalar Ugas)
{

    scalar To = initialTemperatures_[rtiI];

    scalar dTg = Tgas-To;
    const scalar deltaT = mesh_.time().deltaTValue();

    scalar dTe = sqrt(Ugas)/RTI_*(dTg-(1+C_/(sqrt(Ugas)+SMALL))*dTeOld_[rtiI])*deltaT+dTeOld_[rtiI]; 
       
    scalar linkTemperature=initialTemperatures_[rtiI]+dTe;
    dTeOld_[rtiI]=dTe;

    return linkTemperature;
}  

void Foam::rtis::end()
{
    // Do nothing - only valid on write
}


void Foam::rtis::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::rtis::write()
{
    // if (size() && prepare())
    if (size())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }
    writeRTI();
}


Foam::rtis::probeData Foam::rtis::readProbeFile(std::string testFile)
{
    ifstream myfile;
    myfile.open(testFile.c_str(),std::ios::in);
    scalarList time;
    DynamicList<scalarList> data;
    if(myfile.is_open()){
        while(myfile.good()){
            string line;
            getline(myfile,line);
            // skip comment lines
            if(line[0]=='#'){
                continue;
            }
            // skip blank lines
            if(line.size()==0){
                continue;
            }
            std::string word;
            std::istringstream iss(line, std::istringstream::in);
            scalar time1;
            scalar rti1;
            scalarList readRtiValues(size());
            
            iss >> word;
            time1=atof(word.c_str());
            time.append(time1);
            
            label count=0;
            while ( iss >> word ){
                rti1=atof(word.c_str());
                readRtiValues[count]=rti1;
                count++;
            }
            data.append(readRtiValues);
        }

        myfile.close();
    }
    else{
        FatalErrorIn
            (
                "rtis::probeData()"
                )
            << "Unable to open file  \""
            << testFile << "\" for reading. "
            << exit(FatalError);
    }
        
    
    probeData tup(time,data);

    return tup;
}

bool Foam::rtis::interpolateProbeData(const probeData& data)
{
    // see if time exists in probeData
    scalar tBelow=0.0;
    scalar tAbove=0.0;
    label iBelow = -1;
    label iAbove = -1;
    Switch tBelowFound=false;
    Switch tAboveFound=false;
    scalar tSeek = mesh_.time().value();

    for(int i=0;i<data.first().size();i++){
        scalar t=data.first()[i];

        if(t<tSeek){
            tBelow=t;
            iBelow=i;
            tBelowFound=true;
        }
        if(t>=tSeek && !tAboveFound){
            tAbove=t;
            iAbove=i;
            tAboveFound=true;
        }
                
    }

    // if exists, interpolate to find rti probe values for restart
    if(tBelowFound && tAboveFound){
                
        for(int j=0;j<data.second()[iBelow].size();j++){
            scalar valueBelow = data.second()[iBelow][j];
            scalar valueAbove = data.second()[iAbove][j];
            scalar tRatio=(tAbove-tSeek)/(tAbove-tBelow);
            scalar interpolatedValue = -tRatio*(valueAbove-valueBelow)+valueAbove;
            initialTemperatures_[j]=interpolatedValue;
        }

        return true;
    }

    return false;
}
void Foam::rtis::read(const dictionary& dict)
{
    dict.readIfPresent("cleanRestart", cleanRestart_);
    dict.lookup("rtiLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;
    dict.lookup("RTI") >> RTI_;
    dict.lookup("C") >> C_;

    // redetermined all cell locations
    findElements(mesh_);

    if(cleanRestart_){
        dict.lookup("initialTemperatures") >> initialTemperatures_;
        if(initialTemperatures_.size() != (*this).size()){
            FatalErrorIn
                (
                 "rtis::read()"
                 )
                << "Number of initialTemperatures not equal to rtiLocations "
                << "in system/controlDict \""
                << exit(FatalError);
        }
    }
    // get initialTemperatures from RTI probe file
    else{
        initialTemperatures_.resize(size());
        //- locate closest time from <RTI_NAME>/*/RTI file and use those values
        //- loop through time directories searching for match
        fileName rtiSubDir = name_;        
        fileName rtiDir;

        if (Pstream::parRun())
        {
            // Read from undecomposed case
            // (Note: gives problems for distributed data running)
            rtiDir = mesh_.time().path()/".."/rtiSubDir;
        }
        else
        {
            rtiDir = mesh_.time().path()/rtiSubDir;
        }

        // check for directory existence
        struct stat St;
        bool exists = false;
        if (Pstream::master()){
            exists = ( stat(rtiDir.c_str(),&St) == 0 );
        }
        reduce(exists, maxOp<bool>());

        if(!exists){
            WarningIn
                (
                    "rtis::read(const dictionary& dict)"
                    )   << "directory "
                        << rtiDir 
                        << " does not exist"
                        << nl ;
            dict.lookup("initialTemperatures") >> initialTemperatures_;
        }
        else{
            Switch dataNeverFound=true;
            DIR *dirp = opendir(rtiDir.c_str()); 
            
            // Need an array of *dp that can be read in, then sorted
            typedef std::vector<struct dirent *>::iterator direntVectorIterator;
            typedef std::vector<struct dirent *> direntVector;

            direntVector dpVec;

            struct dirent *dp;

            while ((dp = readdir(dirp)) != NULL) 
            {
                dpVec.push_back(dp);
            }

            //  remove "." and ".." directories
            direntVectorIterator iter = dpVec.begin();
            while(iter != dpVec.end()){
                if (strcmp((*iter)->d_name,".")==0 || 
                    strcmp((*iter)->d_name,"..")==0)
                {
                    iter = dpVec.erase(iter);
                }
                else{
                    iter ++;
                }
            }

            //  sort based on time
            std::vector<scalar> time;
            std::vector<scalar>::iterator timeIter;
            iter = dpVec.begin();
            while(iter != dpVec.end()){
                time.push_back(atof((*iter)->d_name));
                iter ++;
            }
            
            // simple bubble sort
            bool sorted=true;

            do{
                sorted=true;
                for(unsigned i=0;i<time.size()-1;i++){
                    if(time[i]>time[i+1]){
                        // swap times and dpVec entries
                        scalar timeTmp = time[i];
                        time[i]=time[i+1];
                        time[i+1]=timeTmp;
                        struct dirent * direntTmp;
                        direntTmp = dpVec[i];
                        dpVec[i]=dpVec[i+1];
                        dpVec[i+1]=direntTmp;
                        sorted=false;
                    }
                }
            }
            while(sorted==false);

            timeIter = time.begin();
            while(timeIter != time.end()){
                timeIter++;
            }

            iter = dpVec.begin();
            while(iter != dpVec.end()){
                            
                fileName testFile=rtiDir/(*iter)->d_name/"RTI";
                
                // read probe data
                probeData data;
                data=readProbeFile(testFile);

                // interpolate probe data to current time
                if(interpolateProbeData(data)){
                    dataNeverFound=false;
                    iter++;
                    break;
                }
                else{
                    iter++;
                    continue;
                }
            }
            (void)closedir(dirp);

            if(dataNeverFound)
            {
                WarningIn
                    (
                        "rtis::read(const dictionary& dict)"
                        )   << "restart data "
                            << " never found"
                            << nl ;
                dict.lookup("initialTemperatures") >> initialTemperatures_;

            }

        }
    }
    
    dict.lookup("activationTemperature") >> activationTemperature_;

    prepare();
}


// ************************************************************************* //
