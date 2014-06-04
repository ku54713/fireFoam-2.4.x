#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "fluidThermoCloud.H"

#include "ThermoSurfaceFilmMeredith.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceFilmModelType(ThermoSurfaceFilmMeredith, fluidThermoCloud);
}
