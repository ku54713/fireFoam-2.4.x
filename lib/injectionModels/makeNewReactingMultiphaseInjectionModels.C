#include "InjectionModel.H"
#include "SprinklerInjection.H"
#include "MultiSprinklerInjection.H"
#include "LookupTableSprinklerInjection.H"

#include "basicReactingMultiphaseCloud.H"

namespace Foam
{
    makeInjectionModelType(SprinklerInjection, basicReactingMultiphaseCloud);
    makeInjectionModelType(MultiSprinklerInjection, basicReactingMultiphaseCloud);
    makeInjectionModelType(LookupTableSprinklerInjection, basicReactingMultiphaseCloud);
}
