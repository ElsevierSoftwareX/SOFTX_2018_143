#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions


# delete old stuff
sh cleanCase.sh
rm -r VTK*

# restore 0
cp -r org 0

runApplication blockMesh

initSurfaceFields

meshUpdater 

foamToVTK -surfaceFields > log.VTK

mv VTK VTK_new

###########################################
# original OF-4-dev dynMesh libs
sh cleanCase.sh
cp -r org 0
runApplication blockMesh 
 
initSurfaceFields 

# ignore doublelinker
echo "
run meshUpdaterOrig"
meshUpdaterOrig  >/dev/null 2>/dev/null |:

echo "
run foamToVTK"
foamToVTK -surfaceFields > log.VTK_Orig

mv VTK VTK_Orig > log.VTK
