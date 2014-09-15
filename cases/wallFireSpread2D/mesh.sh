#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

rm -f log.*

runApplication blockMesh

runApplication setSet -batch system/burner.setSet
runApplication setsToZones -noFlipMap
runApplication createPatch -overwrite

mv -f log.setSet log.setSet.burner
runApplication setSet -batch system/panel.setSet

#runApplication extrudeToRegionMeshNew -AMI -overwrite
runApplication extrudeToRegionMesh -overwrite


# -----------------------------------------------------------------------------
