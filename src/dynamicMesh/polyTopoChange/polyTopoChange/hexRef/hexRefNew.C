/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hexRef.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const bool readHistory = true
)
{
    word hexRefType;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionsD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionsD == 3 && nGeometricD == 2: axisymmetric mesh => hexRefAxi
    // nSolutionsD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionsD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefType = "hexRef";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefType = "hexRefAxi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefType = "hexRef4";
    }

    meshConstructorTable::iterator hexRefIter =
        meshConstructorTablePtr_->find(hexRefType);

    if (hexRefIter == meshContructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionsD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        hexRefIter()(mesh, readHistory)
    );
}

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edge = -GREAT
)
{
    word hexRefType;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionsD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionsD == 3 && nGeometricD == 2: axisymmetric mesh => hexRefAxi
    // nSolutionsD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionsD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefType = "hexRef";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefType = "hexRefAxi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefType = "hexRef4";
    }

    meshConstructorTable::iterator hexRefIter =
        meshConstructorTablePtr_->find(hexRefType);

    if (hexRefIter == meshContructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionsD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        hexRefIter()(mesh, cellLevel, pointLevel, history, level0Edge)
    );
}

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge = -GREAT
)
{
    word hexRefType;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionsD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionsD == 3 && nGeometricD == 2: axisymmetric mesh => hexRefAxi
    // nSolutionsD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionsD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefType = "hexRef";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefType = "hexRefAxi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefType = "hexRef4";
    }

    meshConstructorTable::iterator hexRefIter =
        meshConstructorTablePtr_->find(hexRefType);

    if (hexRefIter == meshContructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionsD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        hexRefIter()(mesh, cellLevel, pointLevel, level0Edge)
    );
}

// ************************************************************************* //
