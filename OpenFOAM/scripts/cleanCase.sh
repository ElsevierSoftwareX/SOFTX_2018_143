#!/bin/sh
# Cleanup case

#remove logFiles
rm -rf log*

# remove all timeSteps 
#rm -rf 0.[0-9]* && rm -rf [1-9]*
rm -rf [0-9]*


# remove Mesh
CASEDIR="$PWD"
MESH="$CASEDIR/constant/polyMesh"
if [ -d $MESH ]
  then
    cd "$MESH"
    rm -rf owner*
    rm -rf neighbour*
    rm -rf point*
    rm -rf face*
    rm -rf boundary*
    rm -rf set*
    rm -rf cell*
    rm -rf level0Edge*
    rm -rf refinementHistory*
    rm -rf surfaceIndex.gz*
    cd $CASEDIR
fi

# remove Processors
rm -rf processor*



