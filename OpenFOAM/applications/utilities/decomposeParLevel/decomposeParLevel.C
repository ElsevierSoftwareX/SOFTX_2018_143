/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.4.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    decomposeParLevel

SourceFiles
    decomposeParLevel.C

Authors
    Christian Kunkelmann and Stefan Batzdorf 2015
    Institute of Technical Thermodynamics
    Technische Universität Darmstadt

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "refinementHistory.H"
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#	include "addTimeOptions.H"
#   include "addRegionOption.H"
#	include "setRootCase.H"
#	include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName regionPrefix = "";
    if (args.options().found("region"))
    {
        regionName = args.options()["region"];
        regionPrefix = regionName;
        Info<< "Operating on region " << regionName << nl << endl;
    }

	polyMesh mesh
	(
		IOobject
		(
			regionName,
			runTime.timeName(),
			runTime
		)
	);

	//- get the number of processors
    int nProcs = 0;
    while
    (
        exists
        (
            args.rootPath()
          / args.caseName()
          / fileName(word("processor") + name(nProcs))
        )
    )
    {
        nProcs++;
    }

	//- read in data for each processor
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        Pout<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(procI))
            << endl;

        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );

        Time& procTime = databases[procI];

        instantList Times = procTime.times();

        //- set startTime and endTime depending on -time and -latestTime options
#       include "checkTimeOptions.H"

        procTime.setTime(Times[startTime], startTime);

        if (procI > 0 && databases[procI-1].value() != procTime.value())
        {
            FatalErrorIn(args.executable())
                << "Time not equal on processors." << nl
                << "Processor:" << procI-1
                << " time:" << databases[procI-1].value() << nl
                << "Processor:" << procI
                << " time:" << procTime.value()
                << exit(FatalError);
        }
    }

    //- check if mesh has already been refined
    //  i.e. if refinementHistory exists

    IOobject historyHeader
    (
        "refinementHistory",
        runTime.timeName(),
        polyMesh::meshSubDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    //- Check if refinementHistory exists
    if (historyHeader.typeHeaderOk<refinementHistory>())
    {
        refinementHistory history (historyHeader);

        labelIOList cellLevel
        (
           IOobject
            (
	    		"cellLevel",
	    		runTime.timeName(),
	    		polyMesh::meshSubDir,
	    		mesh,
	    		IOobject::MUST_READ,
	    		IOobject::NO_WRITE,
	    		false
	    	)
	    );

	    labelIOList pointLevel
	    (
	    	IOobject
	    	(
	    		"pointLevel",
	    		runTime.timeName(),
	    		polyMesh::meshSubDir,
	    		mesh,
	    		IOobject::MUST_READ,
	    		IOobject::NO_WRITE
	    	)
	    );

	    //- distribute information to every processor
        for (label procI = 0; procI < nProcs; procI++)
        {
            labelIOList cellProcAddressing
            (
                IOobject
                (
                    "cellProcAddressing",
                    databases[procI].findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        "cellProcAddressing"
                    ),
                    regionPrefix/polyMesh::meshSubDir,
                    databases[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

	    	labelIOList pointProcAddressing
	    	(
	    		IOobject
	    		(
                    "pointProcAddressing",
                    databases[procI].findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        "pointProcAddressing"
                    ),
                    regionPrefix/polyMesh::meshSubDir,
                    databases[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
	    		)
	    	);

	    	//- update refinementHistory
    		labelList visibleCells (cellProcAddressing.size(), -1);
    		labelList onProc (history.splitCells().size(), 0);

    		forAll(cellProcAddressing, i)
    		{
    			label splitIndex = history.visibleCells()[cellProcAddressing[i]];
    			visibleCells[i] = splitIndex;
    			if (splitIndex != -1)
    			{
    				onProc[splitIndex] = 1;
    				bool foundAllParents = 0;
    				label parentIndex = history.splitCells()[splitIndex].parent_;
    				while (!foundAllParents)
    				{
    					if ((parentIndex != -1) && (onProc[parentIndex] == 0))
    					{
    						onProc[parentIndex] = 1;
    						parentIndex = history.splitCells()[parentIndex].parent_;
    					}
    					else
    					{
    						foundAllParents = 1;
    					}
    				}
    			}
    		}

    		labelList NRemoved (history.splitCells().size(), 0);
    		forAll(NRemoved, i)
    		{
                if (i == 0)
                {
                    if (onProc[0] == 0)
                    {
                        NRemoved[0] = 1;
                    }
                    else
                    {
                        NRemoved[0] = 0;
                    }
                }
                else
                {
                    if (onProc[i] == 0)
        		    {
    	    	        NRemoved[i] = NRemoved[i-1]+1;
    	    	    }
    	    	    else
    	    	    {
    	    	        NRemoved[i] = NRemoved[i-1];
    	    	    }
    	    	}
    		}

    		label NSplitCells = history.splitCells().size();
            if (NRemoved.size() > 0)
            {
                NSplitCells -= NRemoved[NRemoved.size()-1];
            }

    		Info << "Number of split cells overall: " << history.splitCells().size() << " on processor" << procI << ": "
    			 << NSplitCells << " (" << history.splitCells().size() - NSplitCells << " removed)" << endl;

    		forAll (visibleCells, iCell)
    		{
    			if (visibleCells[iCell] != -1)
    			{
    				visibleCells[iCell] -= NRemoved[visibleCells[iCell]];
    			}
    		}

    		refinementHistory::splitCell8 empty;
    		List<refinementHistory::splitCell8> splitCellsProc (NSplitCells, empty);
    		label n = 0;
    		forAll(history.splitCells(), i)
    		{
                if (onProc[i] == 1)
                {
    				splitCellsProc[n] = history.splitCells()[i];

                    if (splitCellsProc[n].parent_ != -1)
                    {
                        splitCellsProc[n].parent_ -= NRemoved[splitCellsProc[n].parent_];
                    }
        			if (history.splitCells()[i].addedCellsPtr_.valid())
    	    		{
                        FixedList<label, 8>& splits = splitCellsProc[n].addedCellsPtr_();
    					forAll(splits, j)
    					{
                            if (splits[j] >= 0)
                            {
                                splits[j] -= NRemoved[splitCellsProc[n].addedCellsPtr_()[j]];
                            }
    					}
    	    		}
    				n++;
                }
    		}

    		//- write refinementHistory
            refinementHistory
            (
                IOobject
                (
                    "refinementHistory",
                    databases[procI].findInstance		//in den selben Ordner schreiben, in dem auch "cellProcAddressing" ist!
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        "cellProcAddressing"
                    ),
                    regionPrefix/polyMesh::meshSubDir,
                    databases[procI],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
    			splitCellsProc,
    			visibleCells,
                true
            ).write();

    		//- write cellLevel
    		labelList cellLevelProc (cellProcAddressing.size(), 0);
    		forAll(cellLevelProc, i)
    		{
    			cellLevelProc[i] = cellLevel[cellProcAddressing[i]];
    		}

            labelIOList
            (
                IOobject
                (
                    "cellLevel",
                    databases[procI].findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        "cellProcAddressing"
                    ),
                    regionPrefix/polyMesh::meshSubDir,
                    databases[procI],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
    			cellLevelProc
            ).write();
	

    		//- write pointLevel
    		labelList pointLevelProc (pointProcAddressing.size(), 0);
    		forAll(pointLevelProc, i)
    		{
    			pointLevelProc[i] = pointLevel[pointProcAddressing[i]];
    		}

            labelIOList
            (
                IOobject
                (
                    "pointLevel",
                    databases[procI].findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        "pointProcAddressing"
                    ),
                    regionPrefix/polyMesh::meshSubDir,
                    databases[procI],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
    			pointLevelProc
            ).write();
        }
    }
    else
    {
        Info << "Refinement history does not exist!" << endl << endl;
    }

	Info << "End" << endl;
    return(0);
}

// ************************************************************************* //
