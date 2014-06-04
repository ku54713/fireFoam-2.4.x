/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "LookupTableVolumeFluxSprinklerInjection.H"
#include "TimeDataEntry.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"
#include "Pstream.H"
#include "OFstream.H"
#include "SortableList.H"

#define DEBUG(x) {                                              \
        std::streamsize p = std::cout.precision();              \
        std::ios::fmtflags myFlags;                             \
        myFlags = cout.flags();                                 \
        std::cout.precision(10);                                \
        std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
        std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
        std::cout << "p" << Pstream::myProcNo();                \
        std::cout << " " << #x " = " << x << std::endl;         \
        std::cout.precision(p);                                 \
        std::cout.flags(myFlags);                               \
    }
#define TRACE(s) {                                              \
        std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
        std::cout << "p" << Pstream::myProcNo();                \
        std::cout << " " << #s << std::endl;                    \
        s;                                                      \
    }

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::LookupTableVolumeFluxSprinklerInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
    )
    :
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    // operatingPressure_(readScalar(this->coeffDict().lookup("operatingPressure"))),
    radiusToSprinkler_(readScalar(this->coeffDict().lookup("radiusToSprinkler"))),
    positionList_(this->coeffDict().lookup("positionList")),
    nSprinklers_(positionList_.size()),
    nActivatedSprinklers_(0),
    injectorCellList_(nSprinklers_),
    tetFaceI_(-1),
    tetPtI_(-1),
    totalParcels_(0),
    direction_(this->coeffDict().lookup("direction")),
    armDirection_(this->coeffDict().lookup("armDirection")),
    parcelDirVec_(direction_),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
        ),
    // flowRateProfile_
    // (
    //     TimeDataEntry<scalar>
    //     (
    //         owner.db().time(),
    //         "flowRateProfile",
    //         this->coeffDict()

    //         )
    //     ),
    tanVec1_(vector::zero),
    tanVec2_(vector::zero),
    Tgas_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            "T"
            )
        ),
    Ugas_
    (
        owner.db().objectRegistry::template lookupObject<volVectorField>
        (
            "U"
            )
        ),
	activeLinks_(this->coeffDict().subDict("rtiCoeffs"). template lookupOrDefault<Switch>("active",false)),
    RTI_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("RTI",200.0)),
    C_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("C",0.0)),
    initialTemperatureList_(nSprinklers_,298.15),
    activationTemperature_(this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("activationTemperature",432.0)),
	// have Andy implement this for Switch in SubModelBase.C
    // activated_(this->template getBaseProperty<Switch>("activated")),
    activatedList_(nSprinklers_,false),
    linkTemperatureList_(nSprinklers_,298.15),
    activationTimeList_(nSprinklers_,GREAT),
    filePtr_(),
    // lookup table related
    parcelIndex_(0),
    activeSprinklerIndex_(-1),
    injectionDeltaT_(0.),
    sampleSize_(0),
    nEle_(0),
    nAzi_(0),
    pressure_(0.),
    kFactor_(0.),
    radius_(0.),
    idealFlowRate_(0.),
    rndGen_(0),
    avgFlux_(),
    dv50_(),
    velMag_(),
    area_(),
    avgVelMag_(),
    ele_(),
    azi_(),
    parcelParticles_(),
    sampleAvgFlux_(),
    sampleDv50_(),
    sampleArea_(),
    sampleAvgVelMag_(),
    sampleEle_(),
    sampleAzi_(),
    sampleParcelParticles_()
{
    
    readTableData();

    idealFlowRate_ = computeIdealFlowRate(); // L/s

    computeInjectionProperties();

    setSampleSizes();

    // scalar injectionDeltaT = 0.001;
    // sampleInjectionTable(injectionDeltaT);
    
    treatSprinklerActivation();

    // Normalise direction vector
    direction_ /= mag(direction_);

    totalParcels_ = sampleSize_;

    // writeVolumeFluxSprinklerInjectionProfile();

    tanVec1_ = armDirection_;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject, gets overwritten later
    // this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_); // m3
    scalar conversionLiterPerM3 = 1000.0;
    this->volumeTotal_ = idealFlowRate_*conversionLiterPerM3*duration_; // m3

    cacheInjectorCells();

    if(activeLinks_){
        writeSprinklerHeader();
    }
}
template<class CloudType>
Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::LookupTableVolumeFluxSprinklerInjection
(
    const LookupTableVolumeFluxSprinklerInjection<CloudType>& im
    )
    :
    InjectionModel<CloudType>(im),
    duration_(im.duration_),
    // operatingPressure_(im.operatingPressure_),
    radiusToSprinkler_(im.radiusToSprinkler_),
    positionList_(im.positionList_),
    nSprinklers_(im.nSprinklers_),
    injectorCellList_(im.injectorCellList_),
    tetFaceI_(im.tetFaceI_),
    tetPtI_(im.tetPtI_),
    totalParcels_(im.totalParcels_),
    direction_(im.direction_),
    armDirection_(im.armDirection_),
    parcelDirVec_(im.parcelDirVec_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    // flowRateProfile_(im.flowRateProfile_),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_),
    Tgas_(im.Tgas_),
    Ugas_(im.Ugas_),
    activeLinks_(im.activeLinks_),
    RTI_(im.RTI_),
    C_(im.C_),
    initialTemperatureList_(im.initialTemperatureList_),
    activationTemperature_(im.activationTemperature_),
	activatedList_(im.activatedList_),
    linkTemperatureList_(im.linkTemperatureList_),
    activationTimeList_(im.activationTimeList_),
    filePtr_(im.filePtr_),
    // lookup table related
    parcelIndex_(im.parcelIndex_),
    activeSprinklerIndex_(im.activeSprinklerIndex_),
    injectionDeltaT_(im.injectionDeltaT_),
    sampleSize_(im.sampleSize_),
    nEle_(im.nEle_),
    nAzi_(im.nAzi_),
    pressure_(im.pressure_),
    kFactor_(im.kFactor_),
    radius_(im.radius_),
    idealFlowRate_(im.idealFlowRate_),
    rndGen_(im.rndGen_),
    avgFlux_(im.avgFlux_),
    dv50_(im.dv50_),
    velMag_(im.velMag_),
    area_(im.area_),
    avgVelMag_(im.avgVelMag_),
    ele_(im.ele_),
    azi_(im.azi_),
    parcelParticles_(im.parcelParticles_),
    sampleAvgFlux_(im.sampleAvgFlux_),
    sampleDv50_(im.sampleDv50_),
    sampleArea_(im.sampleArea_),
    sampleAvgVelMag_(im.sampleAvgVelMag_),
    sampleEle_(im.sampleEle_),
    sampleAzi_(im.sampleAzi_),
    sampleParcelParticles_(im.sampleParcelParticles_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::~LookupTableVolumeFluxSprinklerInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::timeEnd() const
{
	return this->SOI_ + duration_;
}

template<class CloudType>
Foam::scalar Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::timeStart()
{
    nActivatedSprinklers_=0;
    for(label i=0;i<nSprinklers_;i++)
    {
        if(activeLinks_)
        {
            if ( !activatedList_[i])
            {
                computeLinkTemperature(i);
            }
            else
            {
                nActivatedSprinklers_++;
            }
        }
        else
        {
            nActivatedSprinklers_++;
        }
    }

    for(label i=0;i<nSprinklers_;i++)
    {
        std::ostringstream buf;
        buf.precision(2);
        buf << i;

        if
            (
                this->owner().solution().transient()
                && this->owner().db().time().outputTime()
                )
        {
            this->template setBaseProperty<  bool  >
                (word("activated")+buf.str(),activatedList_[i]);
            this->template setBaseProperty< scalar >
                (word("activationTime")+buf.str(),activationTimeList_[i]);
            this->template setBaseProperty< scalar >
                (word("linkTemperature")+buf.str(),linkTemperatureList_[i]);
        }
    }

    if(activeLinks_){
        writeSprinklerData();
    }

    return this->SOI_;
}


template<class CloudType>
Foam::label Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
    )
{

    // set injection time interval for later use
    // injectionDeltaT_ = this->injectionDeltaT_;
    injectionDeltaT_ = time1-time0;
    // DEBUG(injectionDeltaT_);

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        label numberparcels = 
            round((time1 - time0)*(parcelsPerSecond_*nActivatedSprinklers_));
        // DEBUG(parcelsPerSecond_);
        // DEBUG(numberparcels);
        // DEBUG(time1-time0);
        // DEBUG(nActivatedSprinklers_);
        // if(numberparcels == 0)
        //    {return numberparcels;}
        // DEBUG(numberparcels);
        if(numberparcels < totalParcels_*nActivatedSprinklers_)
        {
            return 0;
        }
        else
        {
            // DEBUG(totalParcels_);
            return totalParcels_*nActivatedSprinklers_ ;  // this needs to be modified (kvm)
        }
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
    )
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        if(time1>time0){ 
            // DEBUG("here1");
            // DEBUG(time0);
            // DEBUG(time1);
            // return flowRateProfile_.integrate(time0, time1)*nActivatedSprinklers_;
            return this->volumeTotal_/duration_ * (time1 - time0) * nActivatedSprinklers_;
        }
    }
    else
    {
        return 0.0;
    }
    return 0.0;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setParcelDirVec
(
    scalar elevationAngle, 
    scalar azimuthalAngle
    )
{
    const scalar deg2Rad = pi/180.0;
    scalar alpha = cos(elevationAngle*deg2Rad);
    scalar dCorr = sin(elevationAngle*deg2Rad);
    vector normal = alpha*(tanVec1_*cos(azimuthalAngle*deg2Rad) + tanVec2_*sin(azimuthalAngle*deg2Rad));

    parcelDirVec_ = dCorr*direction_;
    parcelDirVec_ += normal;
    parcelDirVec_ /= mag(parcelDirVec_);
    
    return;
}


/* Called from InjectionModel::inject() in base class */
template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setPositionAndCell
(
    const label parcelI_global,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
    )
{

    // DEBUG("enter setPositionAndCell");

    bool firstParcel=false;
    if(parcelI_global == 0){
        firstParcel = true;
    }

    if(firstParcel){
        // only sample table once per injection
        sampleInjectionTable(injectionDeltaT_);
    }

    // set local parcelI on a per sprinkler basis
    label whichActiveSprinkler=parcelI_global/totalParcels_; 
    // DEBUG(whichActiveSprinkler);
    // If sprinklerIndex = 1, then the second active sprinkler should activate.
    // If sprinklerIndex = 2, then the third active sprinkler should activate.

    activeSprinklerIndex_=-1;
    label activeSprinklerCount=-1;
    for(label si=0;si<activatedList_.size();si++){
        if(activatedList_[si]){
            activeSprinklerCount++;
        }
        if(activeSprinklerCount==whichActiveSprinkler){
            activeSprinklerIndex_=si;
            break;
        }
    }

    label parcelI = parcelI_global%totalParcels_;
    parcelIndex_ = parcelI;

    // ideally, this should be a random value within the min/max of sample angle
    scalar maxEleVar = 2.0;
    scalar maxAziVar = 2.0;
    scalar variationEle = maxEleVar*(rndGen_.scalar01()-0.5);
    scalar variationAzi = maxAziVar*(rndGen_.scalar01()-0.5);
    scalar elevationAngle = sampleEle_[parcelI]+variationEle;
    scalar azimuthalAngle = sampleAzi_[parcelI]+variationAzi;

    // DEBUG(sampleEle_[parcelI]);
    // DEBUG(sampleAzi_[parcelI]);
    
    setParcelDirVec(elevationAngle,azimuthalAngle);

    position = positionList_[activeSprinklerIndex_]+radiusToSprinkler_*parcelDirVec_;
 
    this->findCellAtPosition
        (
         // the first four arguments are returned to InjectionModel::inject()
         cellOwner,
         tetFaceI,
         tetPtI,
         position,
         true
        );

    // DEBUG("exit setPositionAndCell");

    return;
}



template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setParticleDiameter
(
    typename CloudType::parcelType& parcel
)
{

    // DEBUG("enter setParticleDiameter");

    parcel.d() = sampleDv50_[parcelIndex_];

    // DEBUG("exit setParticleDiameter");

    return;
}

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setParticleVelocity
(
    typename CloudType::parcelType& parcel
)
{    
    // DEBUG("enter setParticleVelocity");

    parcel.U() = sampleAvgVelMag_[parcelIndex_] * parcelDirVec_;

    parcel.typeId() = activeSprinklerIndex_; // eventually gets overwritten
    // DEBUG(activeSprinklerIndex_);
    
    // tmp<DimensionedField<scalar, volMesh> > sprinklerId (
    //     new DimensionedField<scalar, volMesh>
    //     (
    //         IOobject
    //         (
    //             this->name() + "sprinklerId" ,
    //             this->owner().time().timeName(),
    //             this->owner(),
    //             IOobject::READ_IF_PRESENT,
    //             IOobject::AUTO_WRITE
    //             ),
    //         this->owner().mesh(),
    //         dimensionedScalar("zero", dimMass, 0.0)
    //         )
    //     );
    // sprinklerId->write();
    
        
    // DEBUG("exit setParticleVelocity");

    return;
}

/* Called from InjectionModel::inject() in base class */
template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
    )
{
    // DEBUG("enter setProperties");

    // set particle diameter
    // DEBUG(parcelI);
    setParticleDiameter(parcel);

    setParticleVelocity(parcel);

    // DEBUG("exit setProperties");
    
    return;
}


template<class CloudType>
bool Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::validInjection(const label)
{
    return true;
}

template<class CloudType>
Foam::scalar Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volume,
    const scalar diameter,
    const scalar rho
    )
{

    scalar nP = sampleParcelParticles_[parcelIndex_];
    
    return nP;
}

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::computeLinkTemperature
(
    const label sprinklerIndex
)
{
 
    scalar To = initialTemperatureList_[sprinklerIndex]; 
    // scalar Tg = 408.0; 
    // scalar U = 1.0;
    // scalar dTg = 110.0;
    // scalar dt = 1.0;
    static List< scalar > dTeOld(nSprinklers_,0.0);
    static List< scalar > linkTemperature(nSprinklers_,To);


    // Initialize values on all procs to -GREAT;
    linkTemperature[sprinklerIndex]=-GREAT;
    scalar U=-GREAT;
    scalar Tg=-GREAT;

    const label ic = injectorCellList_[sprinklerIndex];
    if(ic > -1){
        U=mag(Ugas_[ic])+SMALL;
        Tg=Tgas_[ic];
        scalar dTg = Tg-To;
        const scalar deltaT = this->owner().time().deltaTValue();

        scalar dTe = sqrt(U)/RTI_*(dTg-(1+C_/(sqrt(U)+SMALL))*dTeOld[sprinklerIndex])*deltaT+dTeOld[sprinklerIndex];

        linkTemperature[sprinklerIndex]=initialTemperatureList_[sprinklerIndex]+dTe;
  

        dTeOld[sprinklerIndex]=dTe;
    }
    else{
        dTeOld[sprinklerIndex]=-GREAT;
    }
    // Reduce values to other nodes
    reduce(linkTemperature[sprinklerIndex], maxOp<scalar>());
    reduce(dTeOld[sprinklerIndex], maxOp<scalar>());
    reduce(U, maxOp<scalar>());
    reduce(Tg, maxOp<scalar>());

    linkTemperatureList_[sprinklerIndex] = linkTemperature[sprinklerIndex];
    if(linkTemperature[sprinklerIndex]>activationTemperature_){
        Info << "Sprinkler " << sprinklerIndex << " Activated!\n";
        activatedList_[sprinklerIndex] = true;
        activationTimeList_[sprinklerIndex] = this->owner().time().value();
        this->SOI_ = this->owner().time().value();  
    }

    return;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::cacheInjectorCells()
{
    // Set/cache the injector cell
    // this is the location used for the rti calculation
    forAll(positionList_,index){

        this->findCellAtPosition
            (
                injectorCellList_[index],
                tetFaceI_,
                tetPtI_,
                positionList_[index],
                true
                );

    }
    
    return;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::writeVolumeFluxSprinklerInjectionProfile()
{
    // OFstream osP("sprinklerInjectionProfile");
    
    // osP << "#cell " << "\t";
    // osP << "e1 "    << "\t";
    // osP << "e2 "    << "\t";
    // osP << "a1 "    << "\t";
    // osP << "a2 "    << "\t";
    // osP << "area "  << "\t";
    // osP << "fr "    << "\t";
    // osP << "vfr "   << "\t";
    // osP << "npc "   << "\t";
    // osP << endl;
    
    // for(label iCell = 0; iCell < numberCells; iCell++){
    //     osP << iCell << "\t";
    //     osP << areaEachCell[iCell]  << "\t";
    //     osP << flowRateCell[iCell]*1000 << "\t"; // kg/s
    //     osP << volFCell_[iCell] << "\t";
    //     osP << endl;
    // }

    return;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::treatSprinklerActivation()
{
	if(activeLinks_){
        //sprinkler to be activated via RTI calculation
		this->SOI_=GREAT;
        scalar currentTime = this->owner().time().value();
        
        // read restart properties from <time>/uniform/lagrangian/reactingCloud1/reactingCloud1OutputProperties
        Switch atLeastOneActivated=false;
        for(label i=0;i<nSprinklers_;i++)
        {
            std::ostringstream buf;
            buf.precision(2);
            buf << i;
            
            activatedList_[i] = this->template getBaseProperty<  bool  >
                (word("activated")+buf.str());
            if(activatedList_[i] == true){
                atLeastOneActivated=true;
            }
            activationTimeList_[i] = this->template getBaseProperty< scalar >
                (word("activationTime")+buf.str());
            linkTemperatureList_[i] = this->template getBaseProperty< scalar >
                (word("linkTemperature")+buf.str());
            initialTemperatureList_[i]=linkTemperatureList_[i];
        }

        if(atLeastOneActivated && min(activationTimeList_) < currentTime)
        {
            this->SOI_ = min(activationTimeList_);
        }
        
        for(label i=0;i<nSprinklers_;i++){
            if(initialTemperatureList_[i]==0e0){
                // if linkTemperature not found in reactingCloud1OutputProperties then read from controlDict
                initialTemperatureList_[i] = this->coeffDict().subDict("rtiCoeffs").template lookupOrDefault<scalar>("initialTemperature",298.0);
                activatedList_[i]=false;
                linkTemperatureList_[i] = initialTemperatureList_[i];
                activationTimeList_[i] = GREAT;
            }
        }
	}
    else{
        for(label i=0;i<nSprinklers_;i++){
            activatedList_[i] = true;
        }
    }
 
    return;
}

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::writeSprinklerHeader()
{

    const fileName logDir = "sprinklerPostProcessing"/this->owner().time().timeName();
    mkDir(logDir);
    filePtr_ = new OFstream(logDir/"sprinklerRti.dat");
    OFstream& rtiData = *filePtr_;

    rtiData << "# number of sprinklers: " << nSprinklers_ << endl;
    rtiData << "# positions:" << endl;
    for(label i=0;i<nSprinklers_;i++){
        rtiData << "# sprinkler" << i << " ";
        rtiData << positionList_[i] << endl;
    }
    rtiData << "#t";
    for(label i=0;i<nSprinklers_;i++){
        rtiData << "     ";
        rtiData << "linkTemperature" << i;
    }
    rtiData << endl;

    return;
}

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::writeSprinklerData()
{

    // write to file sprinklerRti.dat
    OFstream& rtiData = *filePtr_;

    rtiData.precision(5);
    rtiData.setf(std::ios::fixed);
    rtiData << this->owner().time().value();
    for(label i=0;i<nSprinklers_;i++){
        rtiData << "     ";
        rtiData << linkTemperatureList_[i];
    }
    rtiData << endl;

    // write to log file
    Info << endl;
    Info << "Sprinkler injection:\n";
    Info << incrIndent;
    Info << indent << "Activated sprinklers: " << nActivatedSprinklers_ << " / " << nSprinklers_ << endl;
    Info << incrIndent;
    for(label i=0;i<nSprinklers_;i++){
        Info << indent << "sprinkler " << i << " (active, link temperature) : " << activatedList_[i] << " " << linkTemperatureList_[i] << endl;
    }
    Info << decrIndent << decrIndent;
    Info << endl;

    return;
}


//- Lookup table functions

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::readTableData()
{

    Info << "Reading sprinkler injection lookup table data\n";

    sampleSize_ = readLabel(this->coeffDict().subDict("lookupTableCoeffs").lookup("sampleSize"));
    tableDirectory_ = this->coeffDict().subDict("lookupTableCoeffs").template lookupOrDefault< word >("tableDirectory","");
    fileName dirName;
    if(Pstream::parRun()){
        dirName = 
            ".."/this->owner().db().time().constant()/tableDirectory_;
    }
    else{
        Info << "tableDirectory_: " << tableDirectory_ << endl;
        dirName = 
            this->owner().db().time().constant()/tableDirectory_;
    }

    {
        IOdictionary dict    
            (        
                IOobject        
                (            
                    "lookup.foam.header",            
                    dirName,
                    // this->owner().db().time().constant()+"/"+tableDirectory_,
                    this->owner().db(),
                    IOobject::MUST_READ_IF_MODIFIED,            
                    IOobject::NO_WRITE        
                    )    
                );
        nEle_ = readLabel(dict.lookup("nEle"));
        nAzi_ = readLabel(dict.lookup("nAzi"));
        pressure_ = readScalar(dict.lookup("pressure"));
        kFactor_ = readScalar(dict.lookup("kFactor"));
        radius_ = readScalar(dict.lookup("radius"));
        Info << "radius " << radius_<< nl;
    }

    {
        IOdictionary dict    
            (        
                IOobject        
                (            
                    "lookup.foam.avgFlux",            
                    dirName,
                    // this->owner().db().time().constant()/tableDirectory_,
                    // this->owner().db().time().constant(),
                    this->owner().db(),
                    IOobject::MUST_READ_IF_MODIFIED,            
                    IOobject::NO_WRITE        
                    )    
                );
        dict.lookup("avgFlux") >> avgFlux_;
    }
     
//    {   
//        IOdictionary dict    
//            (        
//                IOobject        
//                (            
//                    "lookup.foam.dv50",            
//                    dirName,
//                    // this->owner().db().time().constant()/tableDirectory_,
//                    // this->owner().db().time().constant(),
//                    this->owner().db(),
//                    IOobject::MUST_READ_IF_MODIFIED,            
//                    IOobject::NO_WRITE        
//                    )    
//                );
//        dict.lookup("dv50") >> dv50_;
//    }

    {
        IOdictionary dict    
            (        
                IOobject        
                (            
                    "lookup.foam.area",            
                    dirName,
                    // this->owner().db().time().constant()/tableDirectory_,
                    // this->owner().db().time().constant(),
                    this->owner().db(),
                    IOobject::MUST_READ_IF_MODIFIED,            
                    IOobject::NO_WRITE        
                    )    
                );
        dict.lookup("area") >> area_;
    }

//    {
//        IOdictionary dict    
//            (        
//                IOobject        
//                (            
//                    "lookup.foam.avgVelMag",            
//                    dirName,
//                    // this->owner().db().time().constant()/tableDirectory_,
//                    // this->owner().db().time().constant(),
//                    this->owner().db(),
//                    IOobject::MUST_READ_IF_MODIFIED,            
//                    IOobject::NO_WRITE        
//                    )    
//                );
//        dict.lookup("avgVelMag") >> avgVelMag_;
//    }

    {
        IOdictionary dict    
            (        
                IOobject        
                (            
                    "lookup.foam.ele",            
                    dirName,
                    // this->owner().db().time().constant()/tableDirectory_,
                    // this->owner().db().time().constant(),
                    this->owner().db(),
                    IOobject::MUST_READ_IF_MODIFIED,            
                    IOobject::NO_WRITE        
                    )    
                );
        dict.lookup("ele") >> ele_;
    }

    {
        IOdictionary dict    
            (        
                IOobject        
                (            
                    "lookup.foam.azi",            
                    dirName,
                    // this->owner().db().time().constant()/tableDirectory_,
                    // this->owner().db().time().constant(),
                    this->owner().db(),
                    IOobject::MUST_READ_IF_MODIFIED,            
                    IOobject::NO_WRITE        
                    )    
                );
        dict.lookup("azi") >> azi_;
    }

    return;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::computeInjectionProperties()
{

    // compute sprinker orifice diameter
    // K = 29.83 c D^2 
    // where D is in inches, K in british units
    scalar Cd = 0.93; // tapered orifice, orifice coeficient
    scalar conversionKFactor = 14.464; //  gpm/psi^0.5 per lpm/bar^0.5    
    scalar conversionMeterPerInch = 0.0254; //  meter per inch
    scalar factor = 29.83/pow(conversionMeterPerInch,2)*conversionKFactor; // L/min per m2
    Info << "factor: " << factor << endl;
    scalar orificeDiameter = sqrt(kFactor_/factor/Cd); // m
    Info << "orificeDiameter: " << orificeDiameter << endl;

    // compute velocity based on U = mdot/rho/A
    scalar rhow = 1000.0; // density of water, kg/m3
    scalar uJet = kFactor_*sqrt(pressure_)/15.0/rhow/pi/pow(orificeDiameter,2.0); // m
    Info << "uJet: " << uJet << endl;
    velMag_ = uJet*0.8; // account for momentum loss during atomization process

    // compute dv50
    scalar C = 1.9; // sprinkler specific empirical factor for determining dv50 
    scalar sigmaw = 72.8e-3; // N/m @ 20 deg C
    scalar We = rhow*uJet*uJet*orificeDiameter/sigmaw;
    Info << "We: " << We << endl;
    
    dv50_ = C*orificeDiameter/pow(We,0.33333); // m
    Info << "dv50: " << dv50_ << endl;

    // compute particles per parcel
    parcelParticles_.resize(azi_.size());
    forAll(parcelParticles_,nc){
        scalar volumetricFlowRate = avgFlux_[nc]*area_[nc]; // L/s
        scalar injectionDeltaT = 1.0; // s, representative only
        scalar conversionLiterPerM3 = 1000.0;
        scalar volumeToInject = volumetricFlowRate*injectionDeltaT/conversionLiterPerM3; // m3
        // scalar radius = 0.5*dv50_[nc]; // m
        scalar radius = 0.5*dv50_; // m
        scalar singleVolume = 4./3.*pi*pow(radius,3); // m3
        parcelParticles_[nc] = volumeToInject/singleVolume; // None
    }
    
    return;
}


template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::sampleInjectionTable
(
    const scalar injectionDeltaT
)
{

    // scalar specificVolumeToInject = idealFlowRate_*injectionDeltaT;


    // sample twice as many as needed, only keeping top parcelParticles_
    label dSampleSize = 2*sampleSize_;
    List<label> dSampleAziIndicies(dSampleSize);
    List<label> dSampleEleIndicies(dSampleSize);
    SortableList<label> dSampleParcelParticles(dSampleSize);

    forAll(dSampleAziIndicies,i){
        dSampleAziIndicies[i] = rndGen_.integer(0,nAzi_-1);
        dSampleEleIndicies[i] = rndGen_.integer(0,nEle_-1);
        label aziIndex = dSampleAziIndicies[i];
        label eleIndex = dSampleEleIndicies[i];
        label linearIndex = eleIndex+aziIndex*nEle_;
        dSampleParcelParticles[i] = parcelParticles_[linearIndex];
    }
    
    // only keep top half of dSampleParcelParticles
    labelList newToOld;
    sortedOrder(dSampleParcelParticles, newToOld);
    // Info << dSampleParcelParticles << nl;
    // Info << newToOld << nl;

    List<label> sampleAziIndicies(sampleSize_);
    List<label> sampleEleIndicies(sampleSize_);
    // use raw, unsorted values
    if(!false){
        // Set indices to -1, then find random indicie on proc 0, then reduce
        // if other places in the code, like the dispersion model, use random
        // numbers on local processors then the random numbers returned from 
        // rndGen_.integer() will not all be the same here.  Thus, the need
        // for reducing the values to all processors
        forAll(sampleAziIndicies,i){
            sampleAziIndicies[i] = -1;
            sampleEleIndicies[i] = -1;
        }
        if(Pstream::master()){
            forAll(sampleAziIndicies,i){
                sampleAziIndicies[i] = rndGen_.integer(0,nAzi_-1);
                sampleEleIndicies[i] = rndGen_.integer(0,nEle_-1);
            }
        }
        forAll(sampleAziIndicies,i){
            reduce(sampleAziIndicies[i],maxOp<label>());
            reduce(sampleEleIndicies[i],maxOp<label>());
        }
    }
    // use sorted values
    else{
        for(int i=0;i<sampleSize_;i++){
            sampleAziIndicies[i] = dSampleAziIndicies[newToOld[dSampleSize-i-1]];
            sampleEleIndicies[i] = dSampleEleIndicies[newToOld[dSampleSize-i-1]];
            sampleParcelParticles_[i] = dSampleParcelParticles[newToOld[dSampleSize-i-1]];
        }
    }
    // Info << "sampleParcelParticles_" << nl;
    // Info << sampleParcelParticles_ << nl;
    // Info << newToOld << nl;
    
    scalar sumVolumetricFlowRate = 0.0;

    forAll(sampleAziIndicies,nc){
        label aziIndex = sampleAziIndicies[nc];
        label eleIndex = sampleEleIndicies[nc];
        label linearIndex = eleIndex+aziIndex*nEle_;
        sampleAzi_[nc] = azi_[linearIndex];
        sampleEle_[nc] = ele_[linearIndex];
        sampleAvgFlux_[nc] = avgFlux_[linearIndex];
        sampleArea_[nc] = area_[linearIndex];
        // sampleDv50_[nc] = dv50_[linearIndex]; // m
        sampleDv50_[nc] = dv50_; // m
        // sampleAvgVelMag_[nc] = avgVelMag_[linearIndex]; // m/s
        sampleAvgVelMag_[nc] = velMag_; // m/s 

        scalar volumetricFlowRate = sampleAvgFlux_[nc]*sampleArea_[nc];
        sumVolumetricFlowRate += volumetricFlowRate; // L/s
    }

    scalar ratio = sumVolumetricFlowRate/idealFlowRate_;

    // sumVolumetricFlowRate = 0.0;
    // scalar sumParcelParticles = 0.0;

    scalar conversionLiterPerM3 = 1000.0;
    forAll(sampleAziIndicies,nc){
        scalar volumetricFlowRate = sampleAvgFlux_[nc]*sampleArea_[nc]/ratio; // L/s
        // sumVolumetricFlowRate += volumetricFlowRate; // L/s

        scalar volumeToInject = volumetricFlowRate*injectionDeltaT/conversionLiterPerM3; // m3
        scalar radius = 0.5*sampleDv50_[nc]; // m
        scalar singleVolume = 4./3.*pi*pow(radius,3); // m3
        sampleParcelParticles_[nc] = volumeToInject/singleVolume; // None
        // sumParcelParticles += sampleParcelParticles_[nc];
    }
        
    return;
}


template<class CloudType>
Foam::scalar Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::computeIdealFlowRate()
{
    scalar secondsPerMinute = 60.0;
    
    scalar idealFlowRate = kFactor_*sqrt(pressure_); // L/min
    idealFlowRate /= secondsPerMinute; // L/s
    Info << "idealFlowRate " << (idealFlowRate) << " L/s " << nl;

    return idealFlowRate;
}

template<class CloudType>
void Foam::LookupTableVolumeFluxSprinklerInjection<CloudType>::setSampleSizes()
{

    sampleAvgFlux_.resize(sampleSize_);
    sampleDv50_.resize(sampleSize_);
    sampleArea_.resize(sampleSize_);
    sampleAvgVelMag_.resize(sampleSize_);
    sampleEle_.resize(sampleSize_);
    sampleAzi_.resize(sampleSize_);
    sampleParcelParticles_.resize(sampleSize_);
    
    return;
}

// ************************************************************************* //
