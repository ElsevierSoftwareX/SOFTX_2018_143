#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

runApplication rmLog
runApplication rmTimes
runApplication rmMesh
runApplication rmProcessors
runApplication rmPostProcessing

restore0

runApplication blockMesh


# #- Loop to set fields and refine mesh
# maxRefinement=`grep 'maxRefinement' constant/dynamicMeshDict | grep -o -P '(?<=maxRefinement).*(?=;)'`
# for (( c=1; c<=`expr $maxRefinement + 3`; c++ ))
# do

# 	echo "initialize fields and refine iteration: $c"
# 	# initializes an alphaX containing bubble/drop (can create a smooth field) (overwrites timeStep)
# 	initField -x0 0.0 -y0 0.0045 -z0 5e-05 -radius 0.0005 -field alpha.water -latestTime -smooth 2 > log.initField
# #    setFields # would be also possible

# 	a=`expr $maxRefinement + 3`
#     if [[ $c -eq $a ]]; then
#         # if last loop, no further reconstruction necessary, initField maps on the
#         # most recently refined mesh
#         break;
#     fi

# 	# interfaceReconstruct on one TimeStep: use it if in dynamicMesDict a field is used 
# 	# which will be only available after reconstructing the interface e.g. isInterface
# #	reconstructInterface -alphaName alpha.water -latestTime > log.reconstructInterface



# 	# reads all available volFields to refine the mesh as specified in dynamicMeshDict
# 	execRefinement -overwrite #> log.execRefinement 


# done

# rmLog

# rm -r 0/phi

setFields
#- decompose case
#faceMapper

decomposePar 

mpirun -np 2 faceMapper -parallel



touch processor0/case.foam
touch processor1/case.foam



mpirun -np 2 execRefinement -parallel

foamToVTK -case processor0 -surfaceFields > log.VTK
foamToVTK -case processor1 -surfaceFields > log.VTK

#mpirun -np 2 interFoamExtended -parallel > log.

rmLog

# foamToVTK -case processor0 -surfaceFields
# foamToVTK -case processor1 -surfaceFields
# foamToVTK -case processor2 -surfaceFields
# foamToVTK -case processor3 -surfaceFields


# reconstructParMeshAll > log.reconstruct &