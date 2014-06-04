#include "InjectionModel.H"
#include "SprinklerInjection.H"
#include "MultiSprinklerInjection.H"
#include "LookupTableSprinklerInjection.H"

#include "basicThermoCloud.H"

namespace Foam
{
    makeInjectionModelType(SprinklerInjection, basicThermoCloud);
    makeInjectionModelType(MultiSprinklerInjection, basicThermoCloud);
    makeInjectionModelType(LookupTableSprinklerInjection, basicThermoCloud);
}
