#include "InjectionModel.H"
#include "SprinklerInjection.H"
#include "MultiSprinklerInjection.H"
#include "LookupTableSprinklerInjection.H"
#include "LookupTableVolumeFluxSprinklerInjection.H"
#include "UniformSamplingSprinklerInjection.H"

#include "basicReactingCloud.H"

namespace Foam
{
    makeInjectionModelType(SprinklerInjection, basicReactingCloud);
    makeInjectionModelType(MultiSprinklerInjection, basicReactingCloud);
    makeInjectionModelType(LookupTableSprinklerInjection, basicReactingCloud);
    makeInjectionModelType(LookupTableVolumeFluxSprinklerInjection, basicReactingCloud);
    makeInjectionModelType(UniformSamplingSprinklerInjection, basicReactingCloud);
}
