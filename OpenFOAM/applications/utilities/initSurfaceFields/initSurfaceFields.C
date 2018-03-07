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
    initSurfaceFields

Description
    Initializes surfaceFields to test mapping based on face centre position vector

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

    surfaceScalarField mySurfaceScalarField
    (
        IOobject
        (
            "mySurfaceScalarField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mySurfaceScalarField", dimLength, scalar(0))
    );

    surfaceVectorField mySurfaceVectorField
    (
        IOobject
        (
            "mySurfaceVectorField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("mySurfaceVectorField", dimLength, vector(0,0,0))
    );


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

    Info << "Set mySurfaceScalarfield = mag(mesh.Cf())" << endl;
    mySurfaceScalarField = mag(mesh.Cf());

    Info << "Set mySurfaceVectorField = mesh.Cf()" << endl;
    mySurfaceVectorField = mesh.Cf();
    
    runTime.writeNow();

    return(0);
}



// ************************************************************************* //
