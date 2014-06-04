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
  Foam::SprinklerInjection

  Description
  Sprinkler injection

  - User specifies
  - time of start of injection
  - injector position
  - direction (along injection axis)
  - parcel flow rate
  - parcel velocity

  SourceFiles
  SprinklerInjection.C

  \*---------------------------------------------------------------------------*/

#ifndef SprinklerInjection_H
#define SprinklerInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Random.H"
#include "TimeDataEntry.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

    template<class Type>
    class TimeDataEntry;

    class distributionModel;

/*---------------------------------------------------------------------------*\
  Class SprinklerInjection Declaration
  \*---------------------------------------------------------------------------*/

    template<class CloudType>
    class SprinklerInjection
        :
        public InjectionModel<CloudType>
    {
    private:

        // Private data

        //- Injection duration [s]
        const scalar duration_;
        //- Operating pressure [bar]
        const scalar operatingPressure_;
        //spray injection cell deciede by two elevation angles
        scalarList cellEleAngle1_;
        scalarList cellEleAngle2_;
        //spray injection cell deciede by two azimuthal angles
        scalarList cellAziAngle1_;
        scalarList cellAziAngle2_;

        const scalar radiusToSprinkler_;
        // area of each cell 
        // scalarList areaEachCell_;
        // scalarList flowRateCell_;  //lpm
        scalarList nParcelsCell_;
        scalarList volFCell_;


        //- Injector position [m]
        vector position_;

        //- Total number of sprinklers
        label nSprinklers_;

        //- Total number of activated sprinklers
        label nActivatedSprinklers_;

        //- Cell containing injector position []
        label injectorCell_;

        //- Index of tet face for injector cell
        label tetFaceI_;

        //- Index of tet point for injector cell
        label tetPtI_;

        //- Cell containing parcels' position []
        label parcelsCell_;
        //- total number of injected parcels in each time step
        label totalParcels_;
        //- index=1: one parcel at each cell; index=0: parcel number change with flowrate, max 12
        // const label indexParcels_;
        Switch multipleParcelsPerCell_;
        //- Sprinkler Spray Inject Center direction []
        vector direction_;
        //- One Sprinkler Frame Arm direction []
        vector armDirection_;
        vector parcelDirVec_;
// dv50, sigma and gamma for droplet size distribution at each cell
        scalar parcelDv50_;
        scalar parcelSigma_;
        scalar parcelGamma_;
        scalar parcelVelocity_;
        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeDataEntry<scalar> flowRateProfile_;

        //K-factor of operating pressure lpm/bar^0.5
        const label kFactor_;
        //defined range of pressure for fitting coefficients
        scalarList fitPressures_;
        //defined range of speed standard deviation 
        scalarList fitVelocityStdev_;
        //defined azimuthal angles for fitting coefficients
        scalarList fitAziAngles_;
        //fitting coefficients for spray inject flux (lpm/m2) at low pressure
        scalarList fitCoeFluxLowPres_;
        //fitting coefficients for spray inject flux (lpm/m2) at high pressure
        scalarList fitCoeFluxHighPres_;
        scalarList fitCoeFlux_;

        //fitting coefficients for spray volume-median diameter at low pressure
        scalarList fitCoeDv50LowPres_;
        //fitting coefficients for spray volume-median diameter at high pressure
        scalarList fitCoeDv50HighPres_;
        scalarList fitCoeDv50_;
        //fitting coefficients for drop size distribution parameter sigma at low pressure
        scalarList fitCoeSigmaLowPres_;
        //fitting coefficients for drop size distribution parameter sigma at high pressure
        scalarList fitCoeSigmaHighPres_;
        scalarList fitCoeSigma_;
        //fitting coefficients for drop size distribution parameter gamma at low pressure
        scalarList fitCoeGammaLowPres_;
        //fitting coefficients for drop size distribution parameter gamma at high pressure
        scalarList fitCoeGammaHighPres_;
        scalarList fitCoeGamma_;
        //fitting coefficients for drop size-speed at low pressure
        scalarList fitCoeVelocityLowPres_;
        //fitting coefficients for drop size-speed at high pressure
        scalarList fitCoeVelocityHighPres_;
        scalarList fitCoeVelocity_;
        scalarList cellCoeFlux_;
        scalarList cellCoeDv50_;
        scalarList cellCoeSigma_;
        scalarList cellCoeGamma_;
        scalarList cellCoeVelocity_;
        // Tangential vectors to the direction vector

        //- First tangential vector
        vector tanVec1_;

        //- Second tangential vector
        vector tanVec2_;

        const volScalarField& Tgas_;
        const volVectorField& Ugas_;
    	Switch activeLink_;
        scalar RTI_;
        scalar C_;
        scalar initialTemperature_;
        scalar activationTemperature_;
    	Switch activated_;
        label  activatedLabel_;
        scalar linkTemperature_;
        scalar activationTime_;

        //- pointer to log file for sprinkler diagnostic output
        OFstream *filePtr_;

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
        TypeName("sprinklerInjection");


        // Constructors

        //- Construct from dictionary
        SprinklerInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        SprinklerInjection(const SprinklerInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
                (
                    new SprinklerInjection<CloudType>(*this)
                );
        }


        //- Destructor
        virtual ~SprinklerInjection();


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


        void computeLinkTemperature();

        void writeSprinklerHeader();

        void writeSprinklerData();


    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SprinklerInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
