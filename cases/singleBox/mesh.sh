#!/bin/bash -x
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./clean.sh

# create the underlying block mesh
runApplication blockMesh

runApplication surfaceTransformPoints -scale '(42 42 42)' geom/box.stl geom/inches.stl
mv log.surfaceTransformPoints log.surfaceTransformPoints.1
runApplication surfaceTransformPoints -translate '(0 21 0)' geom/inches.stl geom/translated1.stl
mv log.surfaceTransformPoints log.surfaceTransformPoints.2
runApplication surfaceTransformPoints -translate '(0 12 0)' geom/translated1.stl geom/translated2.stl
mv log.surfaceTransformPoints log.surfaceTransformPoints.3
runApplication surfaceTransformPoints -scale '(0.0254 0.0254 0.0254)' geom/translated2.stl geom/meters.stl
mv log.surfaceTransformPoints log.surfaceTransformPoints.4

runApplication snappyHexMesh -overwrite 

runApplication topoSet
mv log.topoSet log.topoSet.box

runApplication topoSet -dict system/topoSetDictBurner
mv log.topoSet log.topoSet.burner

setSet -batch system/createSamplePlane.setSet -time 0 -noZero -constant

runApplication createPatch -overwrite

# Create the wall film region via extrusion
extrudeToRegionMesh -overwrite  -dict system/extrudeToRegionMeshDictFilm 

TMP_FILE=`mktemp /tmp/config.XXXXXXXXXX`
sed -e 's/samplePatch     region0_to_filmRegion/samplePatch     region0_to_pyrolysisRegion/g' constant/filmRegion/polyMesh/boundary > $TMP_FILE
mv $TMP_FILE constant/filmRegion/polyMesh/boundary
sed -e 's/^    region0_to_filmRegion_/    /g' constant/filmRegion/polyMesh/boundary > $TMP_FILE
mv $TMP_FILE constant/filmRegion/polyMesh/boundary

extrudeToRegionMesh -overwrite -dict system/extrudeToRegionMeshDictPyr 


