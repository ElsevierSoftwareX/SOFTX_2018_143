#!/bin/sh

cd ${0%/*} || exit 1



cd src/dynamicMesh/
wmake libso
cd -

cd src/dynamicFvMesh/
wmake libso
cd -

cd applications/utilities/meshUpdater/
wmake
cd - 

cd applications/utilities/initSurfaceFields/
wmake
cd - 
