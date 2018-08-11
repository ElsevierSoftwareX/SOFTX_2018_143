#!/bin/sh
# Author: Daniel Rettenmaier 2018
# Creates a case.foam file for paraview

touch case.foam

if [ -d processor0 ];then
        nproc=`ls -d processor*/ | wc -l`
        for (( i=0; i<=$nproc-1; i++ ))
        do
                touch processor$i/case.foam
        done
fi

