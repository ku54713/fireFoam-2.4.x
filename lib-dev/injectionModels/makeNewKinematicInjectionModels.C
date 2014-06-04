#include "InjectionModel.H"
#include "SprinklerInjection.H"
#include "MultiSprinklerInjection.H"
#include "LookupTableSprinklerInjection.H"

#include "basicKinematicCloud.H"

namespace Foam
{
    makeInjectionModelType(SprinklerInjection, basicKinematicCloud);
    makeInjectionModelType(MultiSprinklerInjection, basicKinematicCloud);
    makeInjectionModelType(LookupTableSprinklerInjection, basicKinematicCloud);
}
