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
    reconstructParLevel

SourceFiles
    reconstructParLevel.C

Authors
    Christian Kunkelmann and Stefan Batzdorf
    Institute of Technical Thermodynamics
    Technische Universität Darmstadt


\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "refinementHistory.H"
#include "DynamicList.H"
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

	labelList cellLevel (mesh.nCells(), 0);
	labelList pointLevel (mesh.nPoints(), 0);
	DynamicList<label> visibleCells;
	DynamicList<refinementHistory::splitCell8> splitCells;

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

    for (label procI = 0; procI < nProcs; procI++)
    {
		//- read refinementHistory
        refinementHistory historyProc
        (
            IOobject
            (
                "refinementHistory",
                databases[procI].findInstance
                (
                    regionPrefix/polyMesh::meshSubDir,
                    "refinementHistory"
                ),
                regionPrefix/polyMesh::meshSubDir,
                databases[procI],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

		label nSplitCells = splitCells.size();

		forAll(historyProc.splitCells(), n)
		{
			refinementHistory::splitCell8 splitCell = historyProc.splitCells()[n];
			if (splitCell.parent_ != -1)
			{
				splitCell.parent_ += nSplitCells;
			}
			if (splitCell.addedCellsPtr_.valid())
			{
				for (int i = 0; i < 8; i++)
				{
					splitCell.addedCellsPtr_()[i] += nSplitCells;
				}
			}
			splitCells.append(splitCell);
		}

		forAll(historyProc.visibleCells(), n)
		{
			label visibleCell = historyProc.visibleCells()[n];
			if (visibleCell != -1)
			{
				visibleCell += nSplitCells;
			}
			visibleCells.append(visibleCell);
		}

		//- read cellProcAddressing
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

		//- read cellLevel
        labelIOList cellLevelProc
        (
            IOobject
            (
                "cellLevel",
                databases[procI].findInstance
                (
                    regionPrefix/polyMesh::meshSubDir,
                    "cellLevel"
                ),
                regionPrefix/polyMesh::meshSubDir,
                databases[procI],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
		
		forAll(cellLevelProc, i)
		{
			cellLevel[cellProcAddressing[i]] = cellLevelProc[i];
		}

		//- read pointProcAddressing
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

		//- read pointLevel
        labelIOList pointLevelProc
        (
            IOobject
            (
                "pointLevel",
                databases[procI].findInstance
                (
                    regionPrefix/polyMesh::meshSubDir,
                    "pointLevel"
                ),
                regionPrefix/polyMesh::meshSubDir,
                databases[procI],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

		forAll(pointLevelProc, i)
		{
			pointLevel[pointProcAddressing[i]] = pointLevelProc[i];
		}

	}

	//- create output
	Info << "creating output files" << endl;

    //- set time
    instantList Times = runTime.times();

    //- set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

	//- write refinementHistory
	List<refinementHistory::splitCell8> tempSplitCells = splitCells;
	List<label> tempVisibleCells (visibleCells.size(), -1);
	forAll (splitCells, iCell)
	{
		tempSplitCells[iCell] = splitCells[iCell];
	}
	forAll (visibleCells, iCell)
	{
		tempVisibleCells[iCell] = visibleCells[iCell];
	}

    refinementHistory
    (
        IOobject
        (
            "refinementHistory",
            runTime.timeName(),
			polyMesh::meshSubDir,
			mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false                       // do not register
        ),
    tempSplitCells,
	tempVisibleCells,
    true
    ).write();

	//- write cellLevel
	labelIOList
	(
		IOobject
		(
			"cellLevel",
			runTime.timeName(),
			polyMesh::meshSubDir,
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
	cellLevel
	).write();

	//- write pointLevel
	labelIOList
	(
		IOobject
		(
			"pointLevel",
			runTime.timeName(),
			polyMesh::meshSubDir,
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
	pointLevel
	).write();

	Info << "End" << endl;

    return(0);
}

// ************************************************************************* //
