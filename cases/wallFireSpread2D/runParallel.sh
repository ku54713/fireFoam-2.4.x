#!/bin/bash

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./clean.sh
./mesh.sh

# Set application name
application="fireFoam"

runApplication decomposePar -force 
mv log.decomposePar log.decomposePar.primaryRegion
runApplication decomposePar -force -region panelRegion
mv log.decomposePar log.decomposePar.panelRegion

#runParallel $application 2
mpirun -np 2 $application -parallel < /dev/null > log.$application 2>log.$application.err

#runApplication reconstructPar


# -----------------------------------------------------------------------------
