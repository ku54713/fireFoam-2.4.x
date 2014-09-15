#!/bin/sh
# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

#runApplication cleanCase
cleanTimeDirectories
runApplication foamClearPolyMesh
mv log.foamClearPolyMesh log.foamClearPolyMesh.primary
runApplication foamClearPolyMesh -region filmRegion
mv log.foamClearPolyMesh log.foamClearPolyMesh.filmRegion
runApplication foamClearPolyMesh -region pyrolysisRegion
mv log.foamClearPolyMesh log.foamClearPolyMesh.pyrolysisRegion

rm -fr *.obj
rm -fr log.*
rm -fr mfr_*
rm -fr cube01_*
rm -fr processor*
rm -fr postProcessing
rm -fr 0/polyMesh
rm -fr 0/filmRegion/cellToRegion
rm -fr 0/pyrolysisRegion/cellToRegion
rm -fr VTK
rm -fr outFlameHeight*
rm -fr *.foam
