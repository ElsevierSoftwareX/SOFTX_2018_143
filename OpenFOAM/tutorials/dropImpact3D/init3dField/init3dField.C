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
    interDyMFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "alphaContactAngleFvPatchScalarField.H"

#include "ReadFields.H"
#include "IOobjectList.H"
#include "DynamicList.H"

#include <iostream>
#include <fstream>

#include "interpolate2dTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // input data
    scalar U0 = 2.0;         // initial droplet veloctiy


#   include "createFields.H"


    forAll(mesh.C(), iCell)
    {
        scalar x = mesh.C()[iCell].component(0);
        scalar y = mesh.C()[iCell].component(1);
        scalar z = mesh.C()[iCell].component(2);
        scalar r = Foam::sqrt(Foam::sqr(x) + Foam::sqr(z));
        scalar phi = Foam::atan(z/(x+VSMALL));

        scalar Urad = interpolate2dTable(r,y,inputUrad,resolution,0.0);
        U[iCell].component(0) = Urad*Foam::cos(phi);
        U[iCell].component(2) = Urad*Foam::sin(phi);

        U[iCell].component(1) = interpolate2dTable(r,y,inputUax,resolution,0.0);

        alpha1[iCell] = interpolate2dTable(r,y,inputAlpha,resolution,0.0);
    }


    // scale (mean) velocity of drop to desired value
    dimensionedVector Udrop = fvc::domainIntegrate(alpha1*U)/fvc::domainIntegrate(alpha1);
    scalar scaleU = -U0/(Udrop.component(1).value()+VSMALL);

    Info << "Scaling U by a factor of " << scaleU << nl << endl;

    U *= scaleU;


    runTime.writeNow();

	Info<< "End\n" << endl;

	return(0);
}


// ************************************************************************* //
