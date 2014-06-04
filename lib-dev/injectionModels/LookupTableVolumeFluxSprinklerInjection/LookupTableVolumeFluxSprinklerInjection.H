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

  Class
  Foam::LookupTableVolumeFluxSprinklerInjection

  Description
  LookupTable-Sprinkler injection

  - User specifies
  - time of start of injection
  - injector positions
  - direction (along injection axis)
  - parcel flow rate
  - parcel velocity

  SourceFiles
  LookupTableVolumeFluxSprinklerInjection.C

  \*---------------------------------------------------------------------------*/

#ifndef LookupTableVolumeFluxSprinklerInjection_H
#define LookupTableVolumeFluxSprinklerInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Random.H"
#include "TimeDataEntry.H"
#include "OFstream.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

    template<class Type>
    class TimeDataEntry;

    class distributionModel;

/*---------------------------------------------------------------------------*\
  Class LookupTableVolumeFluxSprinklerInjection Declaration
  \*---------------------------------------------------------------------------*/

    template<class CloudType>
    class LookupTableVolumeFluxSprinklerInjection
        :
        public InjectionModel<CloudType>
    {
    private:

        // Private data

        //- Injection duration [s]
        const scalar duration_;
        //- Operating pressure [bar]
        // const scalar operatingPressure_;

        const scalar radiusToSprinkler_;

        //- List of injector positions [m]
        List< vector > positionList_;

        //- Total number of sprinklers
        label nSprinklers_;

        //- Total number of activated sprinklers
        label nActivatedSprinklers_;

        //- List of cells containing injector positions []
        List< label > injectorCellList_;

        //- Index of tet face for injector cell
        label tetFaceI_;

        //- Index of tet point for injector cell
        label tetPtI_;

        //- total number of injected parcels in each time step
        label totalParcels_;

        //- Sprinkler Spray Inject Center direction []
        vector direction_;
        //- One Sprinkler Frame Arm direction []
        vector armDirection_;
        vector parcelDirVec_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        // const TimeDataEntry<scalar> flowRateProfile_;

        // Tangential vectors to the direction vector

        //- First tangential vector
        vector tanVec1_;

        //- Second tangential vector
        vector tanVec2_;

        const volScalarField& Tgas_;
        const volVectorField& Ugas_;
        Switch activeLinks_;
        scalar RTI_;
        scalar C_;
        List< scalar > initialTemperatureList_;
        scalar activationTemperature_;
        List< Switch > activatedList_;
        List< label >  activatedLabelList_;
        List< scalar > linkTemperatureList_;
        List< scalar > activationTimeList_;

        //- pointer to log file for sprinkler diagnostic output
        OFstream *filePtr_;


        //- lookup table related
        label parcelIndex_;
        label activeSprinklerIndex_;
        scalar injectionDeltaT_;
        label sampleSize_;
        word tableDirectory_;
        label nEle_;
        label nAzi_;
        scalar pressure_;
        //K-factor of operating pressure lpm/bar^0.5
        scalar kFactor_;
        scalar radius_;
        scalar idealFlowRate_;
        Random rndGen_;

        List< scalar > avgFlux_;
        // List< scalar > dv50_;
        scalar dv50_;
        scalar velMag_;
        List< scalar > area_;
        List< scalar > avgVelMag_;
        List< scalar > ele_;
        List< scalar > azi_;
        List< scalar > parcelParticles_;
        
        List< scalar > sampleAvgFlux_;
        List< scalar > sampleDv50_;
        List< scalar > sampleArea_;
        List< scalar > sampleAvgVelMag_;
        List< scalar > sampleEle_;
        List< scalar > sampleAzi_;
        List< scalar > sampleParcelParticles_;


    protected:

        // sets the diameter of each injected parcel
        void setParticleDiameter(    
            typename CloudType::parcelType& parcel
            );

        // sets the velocity of each injected parcel
        void setParticleVelocity(    
            typename CloudType::parcelType& parcel
            );

        scalar getElevationAngle(Random& rndGen);
        scalar getAzimuthalAngle(Random& rndGen);
        void   setParcelDirVec(scalar elevationAngle,scalar azimuthalAngle);


    public:

        //- Runtime type information
        TypeName("lookupTableVolumeFluxSprinklerInjection");


        // Constructors

        //- Construct from dictionary
        LookupTableVolumeFluxSprinklerInjection
        (
            const dictionary& dict, 
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        LookupTableVolumeFluxSprinklerInjection(const LookupTableVolumeFluxSprinklerInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
                (
                    new LookupTableVolumeFluxSprinklerInjection<CloudType>(*this)
                );
        }


        //- Destructor
        virtual ~LookupTableVolumeFluxSprinklerInjection();


        // Member Functions

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Return the start-of-injection time
        scalar timeStart();

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

        //- Set the injection position and owner cell
        virtual void setPositionAndCell
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            vector& position,
            label& cellOwner,
            label& tetFaceI,
            label& tetPtI
            );

            
        //- Set the parcel properties
        virtual void setProperties
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            typename CloudType::parcelType& parcel
            );

        //- Flag to identify whether model fully describes the parcel
        virtual bool fullyDescribed() const;

        //- Return flag to identify whether or not injection of parcelI is
        //  permitted
        virtual bool validInjection(const label parcelI);


        //- Set number of particles to inject given parcel properties
        virtual scalar setNumberOfParticles
        (
            const label parcels,
            const scalar volume,
            const scalar diameter,
            const scalar rho
            );

        void computeLinkTemperature(const label injectorCell);
        
        void cacheInjectorCells();

        void writeVolumeFluxSprinklerInjectionProfile();

        void treatSprinklerActivation();

        void writeSprinklerHeader();

        void writeSprinklerData();

        //- Table-lookup functions
        void readTableData();
        void computeInjectionProperties();
        void sampleInjectionTable(const scalar injectionDeltaT);
        scalar computeIdealFlowRate();
        void setSampleSizes();


    // Public data

//        //- Phase type enumeration
//        enum lookupDataType
//        {
//            VELOCITY,
//            DIAMETER,
//            VOLUMETRIC_FLUX,
//            UNKNOWN
//        };
//
//        //- Corresponding word representations for lookup data type enumerations
//        static const NamedEnum<lookupDataType, 4> lookupDataTypeNames_;
//
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LookupTableVolumeFluxSprinklerInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
