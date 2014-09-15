#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase
rm -f constant/polyMesh/boundary
rm -f constant/polyMesh/sets/*_old
rm -rf constant/panelRegion/polyMesh
rm -rf patch* 
rm -f log.*
rm -f outFlameHeight*
rm -rf VTK
rm -fr fieldMinMax
rm -fr postProcessing


# -----------------------------------------------------------------------------
