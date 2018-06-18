#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

#clean case
cleanCase.sh

# restore 0
cp -r org/ 0

runApplication blockMesh

# initialize field and set refinement repeatedly (maxRefLevel+1)
for count in {1..3}
do

	echo "initialize fields and refine iteration: $c"
	setFields > log.setFields

	meshUpdater -overwrite
done

# decompose with constraint refinementHistory
runApplication decomposePar
decomposeParLevel

# Make constant/polyMesh folders
for fold in $(ls -d processor*)
do
    mkdir $fold/constant
    cp -r $fold/0/polyMesh $fold/constant
done

# Make .foam file
touch "${PWD##*/}.foam"

#run the solver
echo "run interDyMFoam on 4 procs"
mpirun -np 4 interDyMFoam -parallel > log.interDyMFoam
