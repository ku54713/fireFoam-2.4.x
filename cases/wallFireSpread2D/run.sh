#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./clean.sh
./mesh.sh

# Set application name
application="fireFoam"

$application > log.$application 2> log.error

# -----------------------------------------------------------------------------
