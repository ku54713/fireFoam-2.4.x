#!/bin/bash -x
decomposePar -force >& log.decomposePar.primaryRegion
decomposePar -force -region filmRegion >& log.decomposePar.filmRegion
decomposePar -force -region pyrolysisRegion >& log.decomposePar.pyrolysisRegion
