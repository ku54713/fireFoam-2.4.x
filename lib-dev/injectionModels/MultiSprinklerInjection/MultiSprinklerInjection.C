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

#include "MultiSprinklerInjection.H"
#include "TimeDataEntry.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"
#include "Pstream.H"
#include "OFstream.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiSprinklerInjection<CloudType>::MultiSprinklerInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
    )
    :
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    operatingPressure_(readScalar(this->coeffDict().lookup("operatingPressure"))),
    cellEleAngle1_(this->coeffDict().lookup("cellEleAngle1")),
    cellEleAngle2_(this->coeffDict().lookup("cellEleAngle2")),
    cellAziAngle1_(this->coeffDict().lookup("cellAziAngle1")),
    cellAziAngle2_(this->coeffDict().lookup("cellAziAngle2")),
    radiusToSprinkler_(readScalar(this->coeffDict().lookup("radiusToSprinkler"))),
    nParcelsCell_(cellAziAngle1_.size()),
    volFCell_(cellAziAngle1_.size()),
    positionList_(this->coeffDict().lookup("positionList")),
    nSprinklers_(positionList_.size()),
    nActivatedSprinklers_(0),
    injectorCellList_(nSprinklers_),
    tetFaceI_(-1),
    tetPtI_(-1),
    parcelsCell_(0),
    totalParcels_(0),
    multipleParcelsPerCell_(this->coeffDict().lookup("multipleParcelsPerCell")),
    direction_(this->coeffDict().lookup("direction")),
    armDirection_(this->coeffDict().lookup("armDirection")),
    parcelDirVec_(direction_),
    parcelDv50_(0.000617),
    parcelSigma_(0.784),
    parcelGamma_(0),
    parcelVelocity_(0),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
        ),
    flowRateProfile_
    (
        TimeDataEntry<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()

            )
        ),
    kFactor_
    (
        readScalar(this->coeffDict().lookup("kFactor"))
        ),
    fitPressures_(this->coeffDict().lookup("fitPressureRange")),
    fitVelocityStdev_(this->coeffDict().lookup("fitVelocityStdevRange")),
    fitAziAngles_(this->coeffDict().lookup("fitAzimuthalAngle")),
    fitCoeFluxLowPres_(this->coeffDict().lookup("fitCoeFluxLowPres")),
    fitCoeFluxHighPres_(this->coeffDict().lookup("fitCoeFluxHighPres")),
    fitCoeFlux_(20*(fitAziAngles_.size()-1)+5),
    fitCoeDv50LowPres_(this->coeffDict().lookup("fitCoeDv50LowPres")),
    fitCoeDv50HighPres_(this->coeffDict().lookup("fitCoeDv50HighPres")),
    fitCoeDv50_(20*(fitAziAngles_.size()-1)+5),
    fitCoeSigmaLowPres_(this->coeffDict().lookup("fitCoeSigmaLowPres")),
    fitCoeSigmaHighPres_(this->coeffDict().lookup("fitCoeSigmaHighPres")),
    fitCoeSigma_(20*(fitAziAngles_.size()-1)+5),
    fitCoeGammaLowPres_(this->coeffDict().lookup("fitCoeGammaLowPres")),
    fitCoeGammaHighPres_(this->coeffDict().lookup("fitCoeGammaHighPres")),
    fitCoeGamma_(20*(fitAziAngles_.size()-1)+5),
    fitCoeVelocityLowPres_(this->coeffDict().lookup("fitCoeVelocityLowPres")),
    fitCoeVelocityHighPres_(this->coeffDict().lookup("fitCoeVelocityHighPres")),
    fitCoeVelocity_(20*(fitAziAngles_.size()-1)+5),
    cellCoeFlux_(5*cellEleAngle1_.size()),
    cellCoeDv50_(5*cellEleAngle1_.size()),
    cellCoeSigma_(5*cellEleAngle1_.size()),
    cellCoeGamma_(5*cellEleAngle1_.size()),
    cellCoeVelocity_(5*cellEleAngle1_.size()),
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
    filePtr_()
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
    // Normalise direction vector
    direction_ /= mag(direction_);

    // determine fitting coefficients at operating pressure
    for (label iCell=0; iCell<fitCoeFluxLowPres_.size();iCell++)
    {
        scalar gradCoe;
        gradCoe = (fitCoeFluxHighPres_[iCell]-fitCoeFluxLowPres_[iCell])/(pow(fitPressures_[1],0.5)-pow(fitPressures_[0],0.5));
        fitCoeFlux_[iCell] = fitCoeFluxLowPres_[iCell]+(pow(operatingPressure_,0.5) - pow(fitPressures_[0],0.5))*gradCoe;
        gradCoe = (fitCoeDv50HighPres_[iCell]-fitCoeDv50LowPres_[iCell])/(pow(fitPressures_[1],-1./3.)-pow(fitPressures_[0],-1./3.));
        fitCoeDv50_[iCell] = fitCoeDv50LowPres_[iCell]+(pow(operatingPressure_,-1./3.) - pow(fitPressures_[0],-1./3.))*gradCoe;
        gradCoe = (fitCoeSigmaHighPres_[iCell]-fitCoeSigmaLowPres_[iCell])/(pow(fitPressures_[1],-1./3.)-pow(fitPressures_[0],-1./3.));
        fitCoeSigma_[iCell] = fitCoeSigmaLowPres_[iCell]+(pow(operatingPressure_,-1./3.) - pow(fitPressures_[0],-1./3.))*gradCoe;
        if(kFactor_ == 162){
            gradCoe = (fitCoeGammaHighPres_[iCell]-fitCoeGammaLowPres_[iCell])/(pow(fitPressures_[1],-1./3.)-pow(fitPressures_[0],-1./3.));
            fitCoeGamma_[iCell] = fitCoeGammaLowPres_[iCell]+(pow(operatingPressure_,-1./3.) - pow(fitPressures_[0],-1./3.))*gradCoe;
        }
        gradCoe = (fitCoeVelocityHighPres_[iCell]-fitCoeVelocityLowPres_[iCell])/(pow(fitPressures_[1],0.5)-pow(fitPressures_[0],0.5));
        fitCoeVelocity_[iCell] = fitCoeVelocityLowPres_[iCell]+(pow(operatingPressure_,0.5) - pow(fitPressures_[0],0.5))*gradCoe;
        gradCoe = (fitVelocityStdev_[1]-fitVelocityStdev_[0])/(pow(fitPressures_[1],0.5)-pow(fitPressures_[0],0.5));
        fitVelocityStdev_[0] = fitVelocityStdev_[0]+(pow(operatingPressure_,0.5) - pow(fitPressures_[0],0.5))*gradCoe;
    }

    scalarList fitAllAziAngles(4*(fitAziAngles_.size()-1)+1);

    for(label j=0;j<fitAziAngles_.size();j++){
        fitAllAziAngles[j]=fitAziAngles_[j];
    }
    for (label iQuart=1; iQuart<4;iQuart++)
    {
        for (label jjj=1;jjj<fitAziAngles_.size();jjj++)
        {
            label j = iQuart*(fitAziAngles_.size()-1)+jjj;
            label js = iQuart*(fitAziAngles_.size()-1)-jjj;
            fitAllAziAngles[j]= iQuart*180-fitAllAziAngles[js];
            for (label i=0;i<5;i++)
            {
                label iCellp = j*5+i; label iCelln = js*5+i;
                fitCoeFlux_[iCellp]=fitCoeFlux_[iCelln];
                fitCoeDv50_[iCellp]=fitCoeDv50_[iCelln];
                fitCoeSigma_[iCellp]=fitCoeSigma_[iCelln];
                if(kFactor_ == 162){
                    fitCoeGamma_[iCellp]=fitCoeGamma_[iCelln];
                }
                fitCoeVelocity_[iCellp]=fitCoeVelocity_[iCelln];
            }
        }
    }

    if(!(cellAziAngle1_.size()==cellAziAngle2_.size() && 
         cellAziAngle1_.size()==cellEleAngle1_.size() && 
         cellAziAngle1_.size()==cellEleAngle2_.size()
           )
        )
    {
        Info << nl << "MultiSprinklerInjection cells angle dimension not right" << nl
             << exit(FatalError);
    }
 
    label numberCells = cellAziAngle1_.size();
    scalar totalArea = 0.0;
    scalarList areaEachCell(cellAziAngle1_.size());
    for(label iCell = 0; iCell < numberCells; iCell++)
    {
        scalar theta1 = twoPi*cellEleAngle1_[iCell]/360.0;
        scalar theta2 = twoPi*cellEleAngle2_[iCell]/360.0;
        scalar phi1 = twoPi*cellAziAngle1_[iCell]/360.0;
        scalar phi2 = twoPi*cellAziAngle2_[iCell]/360.0;
        areaEachCell[iCell] = radiusToSprinkler_*radiusToSprinkler_*
            (sin(theta2)-sin(theta1))*(phi2-phi1);
        if(areaEachCell[iCell] < 0.0){
            areaEachCell[iCell] = -areaEachCell[iCell];
        }
        totalArea = totalArea+areaEachCell[iCell];
    }
    for(label i=0;i<5*cellEleAngle1_.size();i++)
    {
        cellCoeFlux_[i]=0.0;
        cellCoeDv50_[i]=0;
        cellCoeSigma_[i]=0;
        cellCoeGamma_[i]=0;
        cellCoeVelocity_[i]=0;
    }

    for(label iCell=0; iCell<numberCells; iCell++)
    {
        scalar iNumCoe(0);
        for (label fitang=0; fitang < (4*(fitAziAngles_.size()-1)+1); fitang ++)
        {
            if(fitAllAziAngles[fitang]>=cellAziAngle1_[iCell] && fitAllAziAngles[fitang]<cellAziAngle2_[iCell])
            {
                iNumCoe++;
                for(label i = 0; i < 5; i++)
                {
                    cellCoeFlux_[i+5*iCell]+=fitCoeFlux_[i+fitang*5];
                    cellCoeDv50_[i+5*iCell]+=fitCoeDv50_[i+fitang*5];
                    cellCoeSigma_[i+5*iCell]+=fitCoeSigma_[i+fitang*5];
                    cellCoeGamma_[i+5*iCell]+=fitCoeGamma_[i+fitang*5];
                    cellCoeVelocity_[i+5*iCell]+=fitCoeVelocity_[i+fitang*5];
                }
            }
        }       
        if(iNumCoe>0){
            for(label i = 0; i < 5; i++)
            {
                cellCoeFlux_[i+5*iCell]/=iNumCoe;
                //Info<<" cell "<<iCell<<" coe flux = "<< cellCoeFlux_[i+5*iCell]<<nl;
                cellCoeDv50_[i+5*iCell]/=iNumCoe;
                cellCoeSigma_[i+5*iCell]/=iNumCoe;
                cellCoeGamma_[i+5*iCell]/=iNumCoe;
                cellCoeVelocity_[i+5*iCell]/=iNumCoe;
            } 
        }
        else{Info<<" No cell coefficients were found !!!!!!!!!!!!!!!!!"<<nl;}
    }    

    scalar totalFlowRate = 0.0; 
    scalarList flowRateCell(cellAziAngle1_.size()); // lpm
    for(label iCell=0; iCell<numberCells; iCell++)
    {
        scalar volumeflux  (0);
        //TODO:  why are eleAngle and j made into integers here?
        scalar eleAngle=int(0.5*(cellEleAngle1_[iCell]+cellEleAngle2_[iCell]));
        scalar aziAngle=int(0.5*(cellAziAngle1_[iCell]+cellAziAngle2_[iCell]));
        if(cellEleAngle2_[iCell] == 90) eleAngle=90;
        if(kFactor_ == 205){
            volumeflux = cellCoeFlux_[5*iCell]+
                cellCoeFlux_[1+5*iCell]*exp(-sqr((eleAngle-15)/7))+
                cellCoeFlux_[2+5*iCell]*exp(-sqr((eleAngle-35)/15))+
                cellCoeFlux_[3+5*iCell]*exp(-sqr((eleAngle-55)/15))+
                cellCoeFlux_[4+5*iCell]*exp(-sqr((eleAngle-90)/10));
        }
        if(kFactor_ == 162){
            volumeflux = cellCoeFlux_[5*iCell]+
                cellCoeFlux_[1+5*iCell]*exp(-sqr((eleAngle-30)/15))+
                cellCoeFlux_[2+5*iCell]*exp(-sqr((eleAngle-45)/15))+
                cellCoeFlux_[3+5*iCell]*exp(-sqr((eleAngle-60)/15))+
                cellCoeFlux_[4+5*iCell]*exp(-sqr((eleAngle-90)/5));
        }
        Info << "cell " << iCell 
             << " a_ang " 
             << aziAngle 
             << " e_ang " << eleAngle;
        Info << " volumeflux " 
             << volumeflux << " flowrate " 
             << areaEachCell[iCell]*volumeflux 
             << nl;
        flowRateCell[iCell] = areaEachCell[iCell]*volumeflux/60/1000.0;  //m^3/s
        if(flowRateCell[iCell] < 0.0){ 
            flowRateCell[iCell] = 0.0;
        }
        totalFlowRate += flowRateCell[iCell]; //m^3/s
    }
    scalar ratioFlowRate = flowRateProfile_.value(1); 
    ratioFlowRate /=totalFlowRate;
    Info << " total flow rate (kg/min) =" << totalFlowRate*1000*60 //kg/min
         << " given rate (kg/min) " << flowRateProfile_.value(1)*1000*60  //kg/min
         << " ratio "<< ratioFlowRate
         << nl;

    scalar minFlowRate = totalFlowRate;
    totalFlowRate = 0.0; 
    for(label iCell=0; iCell<numberCells; iCell++)
    {
        flowRateCell[iCell]*= ratioFlowRate;
        minFlowRate = min(minFlowRate, flowRateCell[iCell]);
        totalFlowRate += flowRateCell[iCell];
    }
    Info << "Correcting flow rates to match theoretical\n";
    ratioFlowRate = flowRateProfile_.value(1); 
    ratioFlowRate /=totalFlowRate;
    Info << " total flow rate (kg/min) =" << totalFlowRate*1000*60 //kg/min
         << " given rate (kg/min) " << flowRateProfile_.value(1)*1000*60  //kg/min
         << " ratio "<< ratioFlowRate
         << nl;

    scalar avgFlowRate =totalFlowRate/(numberCells-1);
    scalar ratioAvgToMin(1);
    ratioAvgToMin= (avgFlowRate-minFlowRate)/6;
    totalParcels_ = 0;
    for(label iCell=0; iCell<numberCells; iCell++)
    {
        nParcelsCell_[iCell]=  int( flowRateCell[iCell]/ratioAvgToMin );
        if(nParcelsCell_[iCell] < 1){
            nParcelsCell_[iCell] = 1;
        }
        if(nParcelsCell_[iCell] > 12){
            nParcelsCell_[iCell] = 12;
        }
        if( ! multipleParcelsPerCell_ ){
            nParcelsCell_[iCell] = 1;
        }
        totalParcels_ += nParcelsCell_[iCell];
        Info << " iCell " << iCell << "  numpcell= " << nParcelsCell_[iCell] << " total= " << totalParcels_ << nl;
    }
    for(label iCell = 0; iCell < numberCells; iCell++)
    {
        volFCell_[iCell] = flowRateCell[iCell]/totalFlowRate ;
    }

    OFstream osP("sprinklerInjectionProfile");
    
    osP << "#cell " << "\t";
    osP << "e1 "    << "\t";
    osP << "e2 "    << "\t";
    osP << "a1 "    << "\t";
    osP << "a2 "    << "\t";
    osP << "area "  << "\t";
    osP << "fr "    << "\t";
    osP << "vfr "   << "\t";
    osP << "npc "   << "\t";
    osP << endl;
    for(label iCell = 0; iCell < numberCells; iCell++){
        osP << iCell << "\t";
        osP << cellEleAngle1_[iCell] << "\t";
        osP << cellEleAngle2_[iCell] << "\t";
        osP << cellAziAngle1_[iCell] << "\t";
        osP << cellAziAngle2_[iCell] << "\t";
        osP << areaEachCell[iCell]  << "\t";
        osP << flowRateCell[iCell]*1000 << "\t"; // kg/s
        osP << volFCell_[iCell] << "\t";
        osP << nParcelsCell_[iCell] << "\t";
        osP << endl;
    }


    tanVec1_ = armDirection_;

    tanVec2_ = direction_^tanVec1_;
    // Set total volume to inject, gets overwritten later
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);

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

    if(activeLinks_){
        writeSprinklerHeader();
    }
}
template<class CloudType>
Foam::MultiSprinklerInjection<CloudType>::MultiSprinklerInjection
(
    const MultiSprinklerInjection<CloudType>& im
    )
    :
    InjectionModel<CloudType>(im),
    duration_(im.duration_),
    operatingPressure_(im.operatingPressure_),
    cellEleAngle1_(im.cellEleAngle1_),
    cellEleAngle2_(im.cellEleAngle2_),
    cellAziAngle1_(im.cellAziAngle1_),
    cellAziAngle2_(im.cellAziAngle2_),
    radiusToSprinkler_(im.radiusToSprinkler_),
    nParcelsCell_(im.nParcelsCell_),
    volFCell_(im.volFCell_),
    positionList_(im.positionList_),
    nSprinklers_(im.nSprinklers_),
    injectorCellList_(im.injectorCellList_),
    tetFaceI_(im.tetFaceI_),
    tetPtI_(im.tetPtI_),
    parcelsCell_(im.parcelsCell_),
    totalParcels_(im.totalParcels_),
    multipleParcelsPerCell_(im.multipleParcelsPerCell_),
    direction_(im.direction_),
    armDirection_(im.armDirection_),
    parcelDirVec_(im.parcelDirVec_),
    parcelDv50_(im.parcelDv50_),
    parcelSigma_(im.parcelSigma_),
    parcelGamma_(im.parcelGamma_),
    parcelVelocity_(im.parcelVelocity_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    kFactor_(im.kFactor_),
    fitPressures_(im.fitPressures_),
    fitVelocityStdev_(im.fitVelocityStdev_),
    fitAziAngles_(im.fitAziAngles_),
    fitCoeFluxLowPres_(im.fitCoeFluxLowPres_),
    fitCoeFluxHighPres_(im.fitCoeFluxHighPres_),
    fitCoeFlux_(im.fitCoeFlux_),
    fitCoeDv50LowPres_(im.fitCoeDv50LowPres_),
    fitCoeDv50HighPres_(im.fitCoeDv50HighPres_),
    fitCoeDv50_(im.fitCoeDv50_),
    fitCoeSigmaLowPres_(im.fitCoeSigmaLowPres_),
    fitCoeSigmaHighPres_(im.fitCoeSigmaHighPres_),
    fitCoeSigma_(im.fitCoeSigma_),
    fitCoeGammaLowPres_(im.fitCoeGammaLowPres_),
    fitCoeGammaHighPres_(im.fitCoeGammaHighPres_),
    fitCoeGamma_(im.fitCoeGamma_),
    fitCoeVelocityLowPres_(im.fitCoeVelocityLowPres_),
    fitCoeVelocityHighPres_(im.fitCoeVelocityHighPres_),
    fitCoeVelocity_(im.fitCoeVelocity_),
    cellCoeFlux_(im.cellCoeFlux_),
    cellCoeDv50_(im.cellCoeDv50_),
    cellCoeSigma_(im.cellCoeSigma_),
    cellCoeGamma_(im.cellCoeGamma_),
    cellCoeVelocity_(im.cellCoeVelocity_),
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
    filePtr_(im.filePtr_)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiSprinklerInjection<CloudType>::~MultiSprinklerInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::timeEnd() const
{
	return this->SOI_ + duration_;
}

template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::timeStart()
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
Foam::label Foam::MultiSprinklerInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
    )
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        label numberparcels = 
            round((time1 - time0)*(parcelsPerSecond_*nActivatedSprinklers_));
        // if(numberparcels == 0)
        //    {return numberparcels;}
        // DEBUG(numberparcels);
        if(numberparcels < totalParcels_*nActivatedSprinklers_)
        {
            return 0;
        }
        else
        {
            return totalParcels_*nActivatedSprinklers_ ;  // this needs to be modified (kvm)
        }
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
    )
{
    // TODO: seems to be a bug here.  should be SOI_ and SOI_+duration_
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        if(time1>time0){ 
            // We never reach here, kvm
            return flowRateProfile_.integrate(time0, time1)*nActivatedSprinklers_;
        }
        else{
            return flowRateProfile_.integrate(time1, time0)*volFCell_[parcelsCell_]/nParcelsCell_[parcelsCell_];
        }
    }
    else
    {
        return 0.0;
    }
}

template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::getElevationAngle(Random& rndGen){
    scalar angle1 = cellEleAngle1_[parcelsCell_];
    scalar angle2 = cellEleAngle2_[parcelsCell_];
    //    scalar i = int(0.5*(angle2 - angle1) + angle1);

    // scalar elevationAngle = rndGen.template sample01<scalar>()*(angle2 - angle1) + angle1;
    scalar randomNumber = rndGen.scalar01(); //bugfix
    // scalar randomNumber = this->owner().rndGen().template sample01<scalar>();
    // randomNumber=0.5;
    // DEBUG(randomNumber);
    scalar elevationAngle = randomNumber*(angle2 - angle1) + angle1;
    //elevationAngle = 0.5*(angle2 - angle1) + angle1;
    return elevationAngle;
}


template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::getAzimuthalAngle(Random& rndGen){
    scalar angle1 = cellAziAngle1_[parcelsCell_];
    scalar angle2 = cellAziAngle2_[parcelsCell_];
    // scalar azimuthalAngle = rndGen.template sample01<scalar>()*(angle2 - angle1) + angle1;
    scalar randomNumber = rndGen.scalar01(); //bugfix
    // randomNumber = this->owner().rndGen().template sample01<scalar>();
    // randomNumber=0.5;
    // DEBUG(randomNumber);
    scalar azimuthalAngle = randomNumber*(angle2 - angle1) + angle1;
    // azimuthalAngle = 0.5*(angle2 - angle1) + angle1;
    //    scalar azimuthalAngle = 0.5*(angle2 - angle1) + angle1;
    return azimuthalAngle;
}


template<class CloudType>
void Foam::MultiSprinklerInjection<CloudType>::setParcelDirVec
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
void Foam::MultiSprinklerInjection<CloudType>::setPositionAndCell
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

    // set local parcelI on a per sprinkler basis
    label whichActiveSprinkler=parcelI_global/totalParcels_; 
    // If sprinklerIndex = 1, then the second active sprinkler should activate.
    // If sprinklerIndex = 2, then the third active sprinkler should activate.
    label activeSprinklerIndex=-1;
    label activeSprinklerCount=-1;
    for(label si=0;si<activatedList_.size();si++){
        if(activatedList_[si]){
            activeSprinklerCount++;
        }
        if(activeSprinklerCount==whichActiveSprinkler){
            activeSprinklerIndex=si;
            break;
        }
    }

    label parcelI = parcelI_global%totalParcels_;

    // the random number generator used to determine the injection location
    //     needs to be used the same number of times on all procs, thus 
    //     maintaining consistency of the injection location across all procs
    // TODO: ensure random number generator for multiSprinklerInjection is 
    // consistent across all processors
    static label seed=0;
    Random rndGen(seed++); // might be more appropriate to use chachedRandom here (static doesn't work for some reason)
    // cachedRandom& rndGen = this->owner().rndGen();

    label numberCells = cellAziAngle1_.size();
    label parcels(0);
    for(label iCell=0;iCell<numberCells;iCell++)
    {
        if(parcelI >= parcels && parcelI < (parcels+nParcelsCell_[iCell]) )
        {
            parcelsCell_=iCell;
        }
        parcels += nParcelsCell_[iCell];
    }

    scalar elevationAngle = getElevationAngle(rndGen);
    scalar azimuthalAngle = getAzimuthalAngle(rndGen);
    
    setParcelDirVec(elevationAngle,azimuthalAngle);

    position = positionList_[activeSprinklerIndex]+radiusToSprinkler_*parcelDirVec_;
 
    this->findCellAtPosition
        (
         // the first four arguments are returned to InjectionModel::inject()
         cellOwner,
         tetFaceI,
         tetPtI,
         position,
         true
        );

    // why are we computing parcelDv50 here in setPositionAndCell?
    if(kFactor_ == 205)
    {
        parcelDv50_ = cellCoeDv50_[5*parcelsCell_]+
            cellCoeDv50_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-15)/7))+
            cellCoeDv50_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-35)/15))+
            cellCoeDv50_[3+5*parcelsCell_]*exp(-sqr((elevationAngle-55)/15))+
            cellCoeDv50_[4+5*parcelsCell_]*exp(-sqr((elevationAngle-90)/10));
        parcelDv50_ *= 0.001;
        parcelSigma_ = cellCoeSigma_[5*parcelsCell_]+
            cellCoeSigma_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-15)/7))+
            cellCoeSigma_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-35)/15))+
            cellCoeSigma_[3+5*parcelsCell_]*exp(-sqr((elevationAngle-55)/15))+
            cellCoeSigma_[4+5*parcelsCell_]*exp(-sqr((elevationAngle-90)/10));
        parcelVelocity_ = cellCoeVelocity_[5*parcelsCell_]+
            cellCoeVelocity_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-15)/10))+
            cellCoeVelocity_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-60)/40));
    }
    if(kFactor_ == 162)
    {
        parcelDv50_ = cellCoeDv50_[5*parcelsCell_]+
            cellCoeDv50_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-30)/15))+
            cellCoeDv50_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-45)/15))+
            cellCoeDv50_[3+5*parcelsCell_]*exp(-sqr((elevationAngle-60)/15))+
            cellCoeDv50_[4+5*parcelsCell_]*exp(-sqr((elevationAngle-90)/5));
        parcelDv50_ *= 0.001; if(parcelDv50_ < 0.0004) parcelDv50_=0.0004;
        parcelSigma_ = cellCoeSigma_[5*parcelsCell_]+
            cellCoeSigma_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-30)/15))+
            cellCoeSigma_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-45)/15))+
            cellCoeSigma_[3+5*parcelsCell_]*exp(-sqr((elevationAngle-60)/15))+
            cellCoeSigma_[4+5*parcelsCell_]*exp(-sqr((elevationAngle-90)/5));
        parcelGamma_ = cellCoeGamma_[5*parcelsCell_]+
            cellCoeGamma_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-30)/15))+
            cellCoeGamma_[2+5*parcelsCell_]*exp(-sqr((elevationAngle-45)/15))+
            cellCoeGamma_[3+5*parcelsCell_]*exp(-sqr((elevationAngle-60)/15))+
            cellCoeGamma_[4+5*parcelsCell_]*exp(-sqr((elevationAngle-90)/5));
        parcelVelocity_ = cellCoeVelocity_[5*parcelsCell_]+
            cellCoeVelocity_[1+5*parcelsCell_]*exp(-sqr((elevationAngle-45)/40));
    }
}



template<class CloudType>
void Foam::MultiSprinklerInjection<CloudType>::setParticleDiameter
(
    typename CloudType::parcelType& parcel
)
{
    label  nBins(61);
    scalarList diameterBins(nBins);
    scalar maxDiameter(4*parcelDv50_);
    maxDiameter =min(0.006, maxDiameter);
    maxDiameter =max(0.001, maxDiameter);
    scalar binWidth=maxDiameter/(nBins-1);

    for(label iBin=0;iBin<nBins;iBin++)
    {
        diameterBins[iBin]=iBin*binWidth;  
    }
    scalarList dropletVolumeFraction(nBins);
    dropletVolumeFraction[0]=0;
    scalar edexp (0);
    for(label iBin=1;iBin<nBins;iBin++)
    {
        scalar d = diameterBins[iBin];
        if(kFactor_ ==205)
        {
            edexp = exp(-sqr(log(d/parcelDv50_))/2/sqr(parcelSigma_)); // Eqn 3-6
            edexp = edexp/(sqrt(twoPi)*parcelSigma_*d)*binWidth;
        }
        if(kFactor_ ==162)
        {
            if(d<=parcelDv50_){
                edexp = exp(-sqr(log(d/parcelDv50_))/2/sqr(parcelSigma_)); // Eqn 3-6
                edexp = edexp/(sqrt(twoPi)*parcelSigma_*d)*binWidth;}
            else {
                edexp = exp(-0.693*pow(diameterBins[iBin-1]/parcelDv50_,parcelGamma_))
                    - exp(-0.693*pow(d/parcelDv50_,parcelGamma_)); }
        }
        dropletVolumeFraction[iBin]=dropletVolumeFraction[iBin-1]+edexp;
    }
    scalar cvfCorrect (1./dropletVolumeFraction[nBins-1]);
    if(cvfCorrect > 1.02)
    {
        for(label iBin=1;iBin<nBins;iBin++)
        {
            scalar d = diameterBins[iBin];
            if(d>parcelDv50_){
                dropletVolumeFraction[iBin]=dropletVolumeFraction[iBin]*cvfCorrect;
            }
        }
    }
    // This rndGen is only executed on the proc that will inject the parcel, 
    // thus causing the rndGen for injection location to be out of sync between procs.
    scalar randomNumber = this->owner().rndGen().template sample01<scalar>();
    // randomNumber=0.5;
    // DEBUG(randomNumber);
    if(randomNumber<0.1){
        randomNumber=0.099;
    }
    if(randomNumber>0.96){
        randomNumber=0.959;
    }
    label iCount = 0;
    while(dropletVolumeFraction[iCount]<randomNumber && iCount<(nBins-1))
    {
        iCount++;  
    }
    if(iCount >= nBins){
        iCount=nBins-1;
    }
    parcel.d()=diameterBins[iCount];

    randomNumber = this->owner().rndGen().template sample01<scalar>();
    // randomNumber=0.5;
    // DEBUG(randomNumber);
    parcel.d() -=  binWidth*randomNumber;

    /* Minimum diameter taken from page 8 of 
       Numerical Simulation of Two Fire Sprinklers Using 
       Measured Spray Characteristics as Starting Conditions Report */
    parcel.d() = max(0.000091,parcel.d());  
    parcel.d() = min(0.006,parcel.d()); 

    return;
}

template<class CloudType>
void Foam::MultiSprinklerInjection<CloudType>::setParticleVelocity
(
    typename CloudType::parcelType& parcel
)
{    
    if(kFactor_ == 205)
    {
        parcelVelocity_ *= (1.-0.8*(exp(-sqr(parcel.d()/0.0005))));
        //    parcelVelocity_ += 2*(this->owner().rndGen().template sample01<scalar>() - 0.5)*fitVelocityStdev_[0];
        scalar randomNumber = this->owner().rndGen().template sample01<scalar>();
        // randomNumber=0.5;
        // DEBUG(randomNumber);
        parcelVelocity_ += randomNumber * fitVelocityStdev_[0];
    }
    if(kFactor_ == 162)
    {
        parcelVelocity_ *= (1.-0.9*(exp(-sqr(parcel.d()/0.0007))));
        //    parcelVelocity_ += 2*(this->owner().rndGen().template sample01<scalar>() - 0.5)*fitVelocityStdev_[0];
        scalar randomNumber = this->owner().rndGen().template sample01<scalar>();
        // randomNumber=0.5;
        // DEBUG(randomNumber);
        parcelVelocity_ += randomNumber * fitVelocityStdev_[0];
    }
    parcel.U() = parcelVelocity_ * parcelDirVec_;
    if(kFactor_ == 162)
    {
        scalar cellAziAngle = 0.5*( cellAziAngle1_[parcelsCell_]+ cellAziAngle2_[parcelsCell_]);
        // These parameters should be inputs coming from the dictionary, not hard coded
        if(cellAziAngle >  -5 && cellAziAngle < 5){
            parcel.d() = 0.0031;
            vector gravityUnitVector = this->owner().g().value()/mag(this->owner().g().value());
            parcel.U() = 2.4*gravityUnitVector;
        }
        if(cellAziAngle > 175 && cellAziAngle < 185){
            parcel.d() = 0.0031;
            vector gravityUnitVector = this->owner().g().value()/mag(this->owner().g().value());
            parcel.U() = 2.4*gravityUnitVector;
        }
    }
    
    return;
}

/* Called from InjectionModel::inject() in base class */
template<class CloudType>
void Foam::MultiSprinklerInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
    )
{
    // set particle diameter
    setParticleDiameter(parcel);

    setParticleVelocity(parcel);

    // Info << "dropletInfo" << " ";
    // Info << parcelI << " ";
    // Info << parcel.d() << " ";
    // Info << parcel.U() << " ";
    // Info << randomNumber << " ";
    // Info << parcel.position() << " ";
    // Info << endl;
    return;
}


template<class CloudType>
bool Foam::MultiSprinklerInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::MultiSprinklerInjection<CloudType>::validInjection(const label)
{
    return true;
}

template<class CloudType>
Foam::scalar Foam::MultiSprinklerInjection<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volume,
    const scalar diameter,
    const scalar rho
    )
{

    // volume passed in here is not valid for sprinkler injection approach
    scalar perParcelVolume = this->volumeToInject(this->injectionDeltaT_,0.0);

    scalar nP = 0.0;
    switch (this->parcelBasis_)
    {
    case MultiSprinklerInjection<CloudType>::pbMass:
    {
        
        // volumeTotal_ only used here
        scalar volume = this->volumeTotal_*nActivatedSprinklers_;
        scalar mass = this->massTotal_*nActivatedSprinklers_;

        nP =
            perParcelVolume/volume
            *mass/rho
            // /(parcels*pi/6.0*pow3(diameter));
            // why in the world does zhoux use this here?!?!?!
            /(pi/6.0*pow3(diameter));
        // Info << "setNumberOfParticles" << " " ;
        // Info << volume << " " ;
        // Info << volumeTotal_ << " " ;
        // Info << massTotal_ << " " ;
        // Info << diameter << " " ;
        // Info << parcels << " " ;
        // Info << nP << " " ;
        // Info << endl;
        break;
    }
    default:
    {
        nP = 0.0;
        FatalErrorIn
            (
                "Foam::scalar "
                "Foam::InjectionModel<CloudType>::setNumberOfParticles"
                "("
                "const label, "
                "const scalar, "
                "const scalar, "
                "const scalar"
                ")"
                )<< "Unknown parcelBasis type" << nl
                 << exit(FatalError);
    }
    }

    return nP;
}

template<class CloudType>
void Foam::MultiSprinklerInjection<CloudType>::computeLinkTemperature
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
void Foam::MultiSprinklerInjection<CloudType>::writeSprinklerHeader()
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
void Foam::MultiSprinklerInjection<CloudType>::writeSprinklerData()
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


// ************************************************************************* //
