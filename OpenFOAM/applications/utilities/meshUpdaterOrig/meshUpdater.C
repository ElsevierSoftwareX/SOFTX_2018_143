/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Application
    meshUpdater

Description
    executes mesh.update() based on constant/dynamicMeshDict
    user has to specify the fields needed for specified refinement criteria in above file

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"

#include "ReadFields.H"
#include "IOobjectList.H"

#include <iostream>
#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // add arguments to function call:
    timeSelector::addOptions(true, false);
#   include "addOverwriteOption.H"
#   include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"


    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    // Read surfaceScalarField (phi)
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    // Read surfaceVectorField (Uf)
    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);


    IOdictionary refineDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar refineInterval = readScalar(refineDict.subDict("dynamicRefineFvMeshCoeffs").lookup("refineInterval"));

    // get current time to reset after mesh update
    scalar time = runTime.time().value();

    for (int i=0; i < refineInterval; i++)
    {
        runTime++;

        mesh.update();
    }

    if(args.optionFound("overwrite"))
    {
        // mesh.update() calls write member if updateTime 
        // if current runTime > first writeTime (of system/controlDict)
        scalar latestTime = runTime.times().last().value();

        // if write was called in mesh.update() delete the new time folder
        if(latestTime > time)
        {
            system("rm -rf "+  runTime.timePath());
        }

        // reset time 
        runTime.setTime(time, 0);
    }

    runTime.writeNow();
    return(0);
}



// ************************************************************************* //
