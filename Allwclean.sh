#!/bin/sh

cd ${0%/*} || exit 1



cd src/dynamicMesh/
wclean
cd -

cd src/dynamicFvMesh/
wclean
cd -

cd applications/utilities/meshUpdaterOrig/
wclean
cd -

cd applications/utilities/meshUpdater/
wclean
cd -

cd applications/utilities/initSurfaceFields/
wclean
cd - 
