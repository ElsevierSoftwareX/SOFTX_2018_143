#!/bin/bash
#
# Author: Daniel Rettenmaier 2018
#
#- Does all the reconstruction steps
#- Spares already reconstructed steps
#- Deletes unfinished reconstructed time steps on abortion with ctrl+c
# 
#- Options: -latestTime, -time 0.01

# checks if a string is contained by a list
containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

# deletes incomplete last time step
function ctrl_c() {
		echo ""
    echo "** Trapped CTRL-C"
    echo ""
  
    if [[ $f0 -eq 1 ]] &&  [[ $f3 -eq 0 ]] ; then
    	rm -r ${timeStep} >/dev/null 2>&1;
    	echo "incomplete time step ${timeStep} deleted"
    fi

    exit
}

#- script options
function options()
{

	local args="$@"

	if [[ $args == *"-latestTime"* ]] ||
	   [[ $args == *"-time"* ]]; then


	        if [[ $args == *"-latestTime"* ]]; then
			# find latestTime for output
			cd processor0
				latestTime=`ls -d [0-9]* | sort -g | tail -1`; 
			cd -
			echo "reconstructing: $latestTime"
		else
			echo "reconstructing: $@"	
		fi

		# deletes -fields entry which is not available nor necessary in reconstructParMesh
		parMeshArg=`echo "$args" | sed 's/-fields[^)]*)//gi'`

		# start reconstruction
		f0=1;f3=0;
		reconstructParMesh $parMeshArg >/dev/null 2>&1; 

    	reconstructPar "$@" >/dev/null 2>&1;

    	reconstructParLevel "$@" >/dev/null 2>&1;
    	f3=1;f0=0;
    	exit
	fi
}

#- reconstruct time single time steps
function reconstructTimeStep()
{
	local tS=$1

	local args="$@"

	# deletes -fields entry which is not available nor necessary in reconstructParMesh
	parMeshArg=`echo "$args" | sed 's/-fields[^)]*)//gi'`

	f0=1;f3=0;
		echo -n " reconstructParMesh $parMeshArg"
        reconstructParMesh $parMeshArg  >/dev/null 2>&1; 

        echo -n " reconstructPar  $args"
    #    reconstructPar "$@" >/dev/null 2>&1 
        grep -q 'ERROR' <<< `reconstructPar "$@" 2>&1` && echo -n "\e[32m ERROR\033[0m"

        echo  " reconstructParLevel  $args"
        reconstructParLevel "$@" >/dev/null 2>&1;
    f3=1;f0=0;
}

#- checks if and which time steps are left for reconstruction
function createTimeListToReconstruct()
{
	echo "Create new List of time steps to reconstruct"
	# find decomposed time steps
	cd processor0
		# get all folders which name begins with a number and sort them
		decompTimeList=`ls -d [0-9]* | sort -g`;
	cd -   >/dev/null 2>&1;

	# find reconstructed time steps
	timeList=`ls -d [0-9]* | sort -g`
	timeList=($timeList)

	#- check if and which time steps are left for reconstruction
	reconstrList=(); timesLeft=false;
	for line in $decompTimeList
	do
		containsElement "$line" "${timeList[@]}"
		reconstruct=$?
		if [ 0 != $reconstruct ]; then
			reconstrList+=("$line");
			timesLeft=true;
		fi
	done
}

#############################################################

SECONDS=0
options $@

#- create new List of time steps
createTimeListToReconstruct

args="$@"
echo "Start: Args: $args"

#- reconstruct time steps and check if new steps are created
while $timesLeft   
do
	count=0;
	for timeStep in ${reconstrList[@]}
	do
		count=$(expr $count + 1);
		echo -e "Reconstruct: $timeStep   \t$count/${#reconstrList[@]}"
		reconstructTimeStep -time $timeStep "${@}"
	done

	sleep 1.5

	#- create new List of time steps
	createTimeListToReconstruct
done


echo "reconstruction finished in "$(($SECONDS / 60)) min and $(($SECONDS % 60)) s .""

# remove temp file
if [ -f logTmp ] ; then
	rm logTmp;
fi




