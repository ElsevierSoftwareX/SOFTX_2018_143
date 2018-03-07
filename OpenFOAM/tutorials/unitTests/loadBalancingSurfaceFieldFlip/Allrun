#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# delete old stuff
sh cleanCase.sh
rm -r VTK
rm -r processor*/VTK*

# restore 0
cp -r org/ 0

# create mesh
runApplication blockMesh

## Try to initialize fields before decomposePar
## it should still flip surfaceVectorFields
# initSurfaceFields

# decompose case
decomposePar

# set surface field in decomposed case
mpirun -np 2 initSurfaceFields -parallel

# do a mesh update, including AMR and LB
mpirun -np 2 meshUpdater -parallel

# preproc for visualization
foamToVTK -surfaceFields -case processor0 -poly > log.VTK
foamToVTK -surfaceFields -case processor1 -poly >> log.VTK

# create case.foam file for paraview on each processor
sh touchCaseFoam.sh
