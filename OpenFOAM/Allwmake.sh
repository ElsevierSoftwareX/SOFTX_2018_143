#!/bin/sh

cd ${0%/*} || exit 1



cd src/dynamicMesh/
wmake libso
cd -

cd src/dynamicFvMesh/
wmake libso
cd -

cd src/twoPhaseProperties/
wmake libso
cd -

cd applications/utilities/meshUpdater/
wmake
cd - 

# linked vs the original OF dynamicMesh library
cd applications/utilities/meshUpdaterOrig/
wmake
cd - 

cd applications/utilities/initSurfaceFields/
wmake
cd - 
