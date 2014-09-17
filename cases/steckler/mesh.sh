#!/bin/bash -x
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./clean.sh

# create the underlying block mesh
runApplication blockMesh

runApplication topoSet -dict system/topoSetDictBurner
mv log.topoSet log.topoSet.burner

runApplication topoSet -dict system/topoSetDictCompartment
mv log.topoSet log.topoSet.compartment

runApplication createPatch -overwrite

# Create 1D and 3D baffles
runApplication createBaffles -overwrite

exit

