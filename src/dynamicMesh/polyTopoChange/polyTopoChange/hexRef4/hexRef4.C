/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "hexRef4.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef4, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<int value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hexRef4::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorIn("hexRef4::reorder(..)") << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


void Foam::hexRef4::getFaceInfo
(
    const label facei,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(facei))
    {
        patchID = mesh_.boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh_.faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


// Adds a face on top of existing facei.
Foam::label Foam::hexRef4::addFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

    label newFacei = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        newFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    return newFacei;
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
Foam::label Foam::hexRef4::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFacei,
    const label meshPointi,
    const face& newFace,
    const label own,
    const label nei
) const
{
    //if (mesh_.isInternalFace(meshFacei))
    //{
    //    return meshMod.setAction
    //    (
    //        polyAddFace
    //        (
    //            newFace,                    // face
    //            own,                        // owner
    //            nei,                        // neighbour
    //            -1,                         // master point
    //            -1,                         // master edge
    //            meshFacei,                  // master face for addition
    //            false,                      // flux flip
    //            -1,                         // patch for face
    //            -1,                         // zone for face
    //            false                       // face zone flip
    //        )
    //    );
    //}
    //else
    {
        // Two choices:
        // - append (i.e. create out of nothing - will not be mapped)
        //   problem: field does not get mapped.
        // - inflate from point.
        //   problem: does interpolative mapping which constructs full
        //   volPointinterpolation!

        // For now create out of nothing

        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );


        ////- Inflate-from-point:
        //// Check if point has any internal faces we can use.
        //label masterPointi = -1;
        //
        //const labelList& pFaces = mesh_.pointFaces()[meshPointi];
        //
        //forAll(pFaces, i)
        //{
        //    if (mesh_.isInternalFace(pFaces[i]))
        //    {
        //        // meshPoint uses internal faces so ok to inflate from it
        //        masterPointi = meshPointi;
        //
        //        break;
        //    }
        //}
        //
        //return meshMod.setAction
        //(
        //    polyAddFace
        //    (
        //        newFace,                    // face
        //        own,                        // owner
        //        nei,                        // neighbour
        //        masterPointi,               // master point
        //        -1,                         // master edge
        //        -1,                         // master face for addition
        //        false,                      // flux flip
        //        -1,                         // patch for face
        //        -1,                         // zone for face
        //        false                       // face zone flip
        //    )
        //);
    }
}


// Modifies existing facei for either new owner/neighbour or new face points.
void Foam::hexRef4::modFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[facei])
     || (
            mesh_.isInternalFace(facei)
         && (nei != mesh_.faceNeighbour()[facei])
        )
     || (newFace != mesh_.faces()[facei])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    facei,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    facei,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


// Bit complex way to determine the unrefined edge length.
Foam::scalar Foam::hexRef4::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::getLevel0EdgeLength() const"
        )   << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, celli)
        {
            const label cLevel = cellLevel_[celli];

            const labelList& cEdges = mesh_.cellEdges(celli);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    forAll(cellLevel_, celli)
    {
        const label cLevel = cellLevel_[celli];

        const labelList& cEdges = mesh_.cellEdges(celli);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == GREAT2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexRef4::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorIn("hexRef4::getLevel0EdgeLength()")
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


// Check whether pointi is an anchor on celli.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef4::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    if (cellAnchorPoints[celli].size())
    {
        label index = findIndex(cellAnchorPoints[celli], pointi);

        if (index != -1)
        {
			if (index >= 4) //AB....
			{
				if (index == 4)
				{
					index = 8;
				}
				index = 8 - index;
			} //AB
            return cellAddedCells[celli][index];
        }


        // pointi is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[facei];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[celli], f[fp]);

            if (index != -1)
            {
                if (index >= 4) //AB....
                {
                    if (index == 4)
                    {
                        index = 8;
                    }
                    index = 8 - index;
                } //AB
                return cellAddedCells[celli][index];
            }
        }

        // Problem.
        dumpCell(celli);
        Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
            << endl;

        FatalErrorIn("hexRef4::getAnchorCell(..)")
            << "Could not find point " << pointi
            << " in the anchorPoints for cell " << celli << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return celli;
    }
}


// Get new owner and neighbour
void Foam::hexRef4::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label facei,
    const label pointi,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[facei],
        facei,
        pointi
    );

    if (mesh_.isInternalFace(facei))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[facei],
            facei,
            pointi
        );
    }
    else
    {
        nei = -1;
    }
}


// Get point with the lowest pointLevel
Foam::label Foam::hexRef4::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


// Get point with the highest pointLevel
Foam::label Foam::hexRef4::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexRef4::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexRef4::dumpCell(const label celli) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(celli) + ".obj");
    Pout<< "hexRef4 : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[celli];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointi = f[fp];
            label nexPointi = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointi]+1
                << ' ' << pointToObjVert[nexPointi]+1 << nl;
        }
    }
}


// Find point with certain pointLevel. Skip any higher levels.
Foam::label Foam::hexRef4::findLevel
(
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointi = f[fp];

        if (pointLevel_[pointi] < wantedLevel)
        {
            FatalErrorIn("hexRef4::findLevel(..)")
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointi] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    FatalErrorIn("hexRef4::findLevel(..)")
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


// Gets cell level such that the face has four points <= level.
Foam::label Foam::hexRef4::getAnchorLevel(const label facei) const
{
    const face& f = mesh_.faces()[facei];

    if (f.size() <= 4)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[facei]];

        if (countAnchors(f, ownLevel) == 4)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 4)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::hexRef4::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(neiPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&n) / (dir&n);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexRef4::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(boundaryPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkBoundaryOrientation(..)")
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningIn("checkBoundaryOrientation(..)")
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
            //<< abort(FatalError);
    }
}


// If p0 and p1 are existing vertices check if edge is split and insert
// splitPoint.
void Foam::hexRef4::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


// Internal faces are one per edge between anchor points. So one per midPoint
// between the anchor points. Here we store the information on the midPoint
// and if we have enough information:
// - two anchors
// - two face mid points
// we add the face. Note that this routine can get called anywhere from
// two times (two unrefined faces) to four times (two refined faces) so
// the first call that adds the information creates the face.
Foam::label Foam::hexRef4::storeMidPointinfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const bool faceOrder,
    const label edgeMidPointi,
    const label anchorPointi,
    const label faceMidPointi,

    Map<edge>& midPointToAnchors,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointi);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointi, edge(anchorPointi, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointi != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointi;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }


    // Check if this call of storeMidPointinfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors)
    {
        const cell& cFaces = mesh_.cells()[celli];
        label face1 = -1;

        forAll(cFaces, i)
        {
            label faceJ = cFaces[i];

            if (faceMidPoint[faceJ] != faceMidPointi
                && faceMidPoint[faceJ] >= 0
                && faceMidPoint[faceJ] != 1234567890) // could be replace by faceJ == type "empty"
            {
                face1 = faceJ;  // needed to take the other face points
            }
        }

        const edge& anchors = midPointToAnchors[edgeMidPointi];
        label index = findIndex(cellAnchorPoints[celli], anchorPointi);

        if (findIndex(cellAnchorPoints[celli], anchorPointi) == 0)
        {
            index = 4;
        }
        if (findIndex(cellAnchorPoints[celli], anchorPointi) == 4)
        {
            index = 8;
        }

        label point1 = cellAnchorPoints[celli][8 - index];
        label edgeMidPointJ = -1;

        const face& f = mesh_.faces()[face1];       //other face
        const labelList& fEdges = mesh_.faceEdges(face1);

        DynamicList<label> newFaceVerts(4);

        if (faceOrder == (mesh_.faceOwner()[facei] == celli))
        {
            label anch = findIndex(f, point1);

            if (pointLevel_[f[f.rcIndex(anch)]] <= cellLevel_[celli])
            {
                label edgeJ = fEdges[f.rcIndex(anch)];
                edgeMidPointJ = edgeMidPoint[edgeJ];
            }
            else
            {
                label edgeMid = findLevel(f, f.rcIndex(anch), false, cellLevel_[celli] +1);
                edgeMidPointJ = f[edgeMid];
            }

            newFaceVerts.append(faceMidPointi);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointi,  // edge between faceMid and
                edgeMidPointi,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                edgeMidPointJ,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointJ);
            newFaceVerts.append(faceMidPoint[face1]);
        }
        else
        {
            label anch = findIndex(f, point1);

            if (pointLevel_[f[f.fcIndex(anch)]] <= cellLevel_[celli])
            {
                label edgeJ = fEdges[anch];
                edgeMidPointJ = edgeMidPoint[edgeJ];
            }
            else
            {
                label edgeMid = findLevel(f, f.fcIndex(anch), false, cellLevel_[celli] +1);
                edgeMidPointJ = f[edgeMid];
            }

            newFaceVerts.append(edgeMidPointJ);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointJ,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                faceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointi);
            newFaceVerts.append(faceMidPoint[face1]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchorPointi
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchors.otherVertex(anchorPointi)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointi];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            neiPt = mesh_.points()[anchorPointi];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointi];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                neiPt = mesh_.points()[anchorPointi];
            }

            checkInternalOrientation
            (
                meshMod,
                celli,
                facei,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            facei,
            anchorPointi,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


// Creates all the 4 internal faces for celli.
void Foam::hexRef4::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[celli];
    const label cLevel = cellLevel_[celli];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(8);
    // From edge mid to face mids
    Map<edge> midPointToFaceMids(8);

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        const face& f = mesh_.faces()[facei];
        const labelList& fEdges = mesh_.faceEdges(facei, storage);

        // We are on the celli side of face f. The face will have 1 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointi = -1;

        label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.
            // In 2D case this should never happen, as the shared faces with uncutted mesh
            // can be only the ones splitted in two which have two anchors and two middle edge points.

            Info << "Should never happen: nAnchors == 1" << endl;

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel+1 point: 
            // these two lines just repeat the same steps two times to reach the second cLevel+1
            label edgeMid = findLevel(f, f.fcIndex(anchorFp), true, cLevel+1);
            // to be checked!!!
            label faceMid = findLevel(f, f.fcIndex(edgeMid), true, cLevel+1);

            faceMidPointi = f[faceMid];
        }
        else if (nAnchors == 2)
        {
            // Only two anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            faceMidPointi = 1234567890;        // again not nice

        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting. For non-empty faces, faceMidPoint[facei] will be 1234567890,
            // which is not a point ID. This is checked in storeMidPointinfo.

            faceMidPointi = faceMidPoint[facei];
        }
        else
        {
            FatalErrorIn("createInternalFaces")
                << "nAnchors:" << nAnchors
                << " facei:" << facei
                << abort(FatalError);
        }



        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointinfo will collect
        // all the info and combine it all.
        // Only empty faces are considered as they are split like the 3D case.

        if(faceMidPoint[facei] != 1234567890 && faceMidPointi != 1234567890)
        {
            forAll(f, fp0)
            {
                label point0 = f[fp0];

                if (pointLevel_[point0] <= cLevel)
                {
                    // Anchor.

                    // Walk forward
                    // ~~~~~~~~~~~~
                    // to cLevel+1 or edgeMidPoint of this level.


                    label edgeMidPointi = -1;

                    label fp1 = f.fcIndex(fp0);

                    if (pointLevel_[f[fp1]] <= cLevel)
                    {
                        // Anchor. Edge will be split.
                        label edgeI = fEdges[fp0];

                        edgeMidPointi = edgeMidPoint[edgeI];

                        if (edgeMidPointi == -1)
                        {
                            dumpCell(celli);

                            const labelList& cPoints = mesh_.cellPoints(celli);

                            FatalErrorIn("createInternalFaces(..)")
                                << "cell:" << celli << " cLevel:" << cLevel
                                << " cell points:" << cPoints
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, cPoints)()
                                << " face:" << facei
                                << " f:" << f
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, f)()
                                << " faceAnchorLevel:" << faceAnchorLevel[facei]
                                << " faceMidPoint:" << faceMidPoint[facei]
                                << " faceMidPointi:" << faceMidPointi
                                << " fp:" << fp0
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // Search forward in face to clevel+1
                        // In 2D this should never be used as this case is verified with
                        // shared faces between refined and not-refined cells which have 
                        // to be refined in this turn.

                        label edgeMid = findLevel(f, fp1, true, cLevel+1);

                        edgeMidPointi = f[edgeMid];
                    }

                    label newFacei = storeMidPointinfo
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        cellMidPoint,
                        faceMidPoint,
                        edgeMidPoint,

                        celli,
                        facei,
                        true,                   // mid point after anchor
                        edgeMidPointi,          // edgemid
                        point0,                 // anchor
                        faceMidPointi,

                        midPointToAnchors,
                        meshMod
                    );

                    if (newFacei != -1)
                    {
                        nFacesAdded++;

                        if (nFacesAdded == 4)
                        {
                            break;
                        }
                    }



                    // Walk backward
                    // ~~~~~~~~~~~~~

                    label fpMin1 = f.rcIndex(fp0);

                    if (pointLevel_[f[fpMin1]] <= cLevel)
                    {
                        // Anchor. Edge will be split.
                        label edgeI = fEdges[fpMin1];

                        edgeMidPointi = edgeMidPoint[edgeI];

                        if (edgeMidPointi == -1)
                        {
                            dumpCell(celli);

                            const labelList& cPoints = mesh_.cellPoints(celli);

                            FatalErrorIn("createInternalFaces(..)")
                                << "cell:" << celli << " cLevel:" << cLevel
                                << " cell points:" << cPoints
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, cPoints)()
                                << " face:" << facei
                                << " f:" << f
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, f)()
                                << " faceAnchorLevel:" << faceAnchorLevel[facei]
                                << " faceMidPoint:" << faceMidPoint[facei]
                                << " faceMidPointi:" << faceMidPointi
                                << " fp:" << fp0
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // Search back to clevel+1
                        // In 2D this should never be used as this case is verified with
                        // shared faces between refined and not-refined cells which have 
                        // to be refined in this turn.

                        label edgeMid = findLevel(f, fpMin1, false, cLevel+1);

                        edgeMidPointi = f[edgeMid];
                    }

                    newFacei = storeMidPointinfo
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        cellMidPoint,
                        faceMidPoint,
                        edgeMidPoint,

                        celli,
                        facei,
                        false,                  // mid point before anchor
                        edgeMidPointi,          // edgemid
                        point0,                 // anchor
                        faceMidPointi,

                        midPointToAnchors,
                        meshMod
                    );

                    if (newFacei != -1)
                    {
                        nFacesAdded++;

                        if (nFacesAdded == 4)
                        {
                            break;
                        }
                    }
                }   // done anchor
            } // done face
        }   // done check empty condition

        if (nFacesAdded == 4)
        {
            break;
        }
    }
}


void Foam::hexRef4::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


// Same as walkFaceToMid but now walk back.
void Foam::hexRef4::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::hexRef4::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


// Debug: check if wanted refinement is compatible with 2:1
void Foam::hexRef4::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorIn
            (
                "hexRef4::checkWantedRefinementLevels(const labelList&)"
            )   << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label facei = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchI = mesh_.boundaryMesh().whichPatch(facei);

            dumpCell(own);
            FatalErrorIn
            (
                "hexRef4::checkWantedRefinementLevels(const labelList&)"
            )   << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << facei
                << " on patch " << patchI << " "
                << mesh_.boundaryMesh()[patchI].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


// Set instance for mesh files
void Foam::hexRef4::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexRef4::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef4::hexRef4(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("level0Edge", dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (readHistory ? mesh_.nCells() : 0)  // All cells visible if not be read
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&)"
        )   << "History enabled but number of visible cells "
            << history_.visibleCells().size() << " in "
            << history_.objectPath()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&)"
        )   << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "cellLevel file " << cellLevel_.objectPath() << endl
            << "pointLevel file " << pointLevel_.objectPath() << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        history
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&, const refinementHistory&)"
        )   << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&, const refinementHistory&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        List<refinementHistory::splitCell8>(0),
        labelList(0),
        false
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Get points of a cell (without using cellPoints addressing)
Foam::labelList Foam::hexRef4::cellPoints(const label celli) const
{
    // Pick up points of the cell
    const cell& cFaces = mesh_.cells()[celli];

    labelHashSet cPoints(4*cFaces.size());

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            cPoints.insert(f[fp]);
        }
    }
    return cPoints.toc();
}


Foam::labelList Foam::hexRef4::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


// Given a list of cells to refine determine additional cells to refine
// such that the overall refinement:
// - satisfies maxFaceDiff (e.g. 2:1) across neighbouring faces
// - satisfies maxPointDiff (e.g. 4:1) across selected point connected
//   cells. This is used to ensure that e.g. cells on the surface are not
//   point connected to cells which are 8 times smaller.
Foam::labelList Foam::hexRef4::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
)
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef4::consistentSlowRefinement"
            "(const label, const labelList&, const labelList&"
            ", const label, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialize it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, celli)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[celli] = refinementData
        (
            maxFaceDiff*(cellLevel_[celli]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[celli]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].count() = allCellInfo[celli].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef4::consistentSlowRefinement"
                "(const label, const labelList&, const labelList&"
                ", const label, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[facei]];

        if (mesh_.isInternalFace(facei))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[facei]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces inbetween different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];
            label nei = faceNeighbour[facei];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                facei,
                allCellInfo[own],
                FaceCellWave<refinementData>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(facei);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexRef4::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointi)
        {
            label& pLevel = maxPointCount[pointi];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                refinementData& cellInfo = allCellInfo[celli];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointi]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            facei,
                            celli,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[facei].count())
                        {
                            changedFacesInfo.insert(facei, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorIn
                (
                    "hexRef4::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << facei << " faceRefData:" << allFaceInfo[facei]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchI = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorIn
                (
                    "hexRef4::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << facei
                    << " refData:" << allFaceInfo[facei]
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef4::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef4::consistentSlowRefinement2"
            "(const label, const labelList&, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[celli],
            cellLevel_[celli]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, celli)
    {
        if (!allCellInfo[celli].valid(dummyTrackData))
        {
            allCellInfo[celli] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[celli],
                cellLevel_[celli]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef4::consistentSlowRefinement2"
                "(const label, const labelList&, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }

        label own = faceOwner[facei];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(facei))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[facei];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[facei].updateFace
            (
                mesh_,
                facei,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(facei);
        seedFacesInfo.append(allFaceInfo[facei]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces inbetween different refinement levels for now
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    //if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().timeName(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar("zero", dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, celli)
    //    {
    //        wantedLevel[celli] = allCellInfo[celli].wantedLevel(cc[celli]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].originLevel() = sizeof(label)*8-2;
        allCellInfo[celli].origin() = cc[celli];
    }

    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, celli)
    {
        label wanted = allCellInfo[celli].wantedLevel(cc[celli]);

        if (wanted > cellLevel_[celli]+1)
        {
            refineCell.set(celli);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
//        if (refineCell.get(celli))
        if (refineCell[celli])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
//        if (refineCell.get(celli))
        if (refineCell[celli])
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli))
                {
                    cellsOut2.insert(celli);
                }
            }
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli) && !savedRefineCell.get(celli))
                {
                    dumpCell(celli);
                    FatalErrorIn
                    (
                        "hexRef4::consistentSlowRefinement2"
                        "(const label, const labelList&, const labelList&)"
                    )   << "Cell:" << celli << " cc:"
                        << mesh_.cellCentres()[celli]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}


// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRef4::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();


    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, celli)
    {
        newCellLevel.append(cellLevel_[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointi)
    {
        newPointLevel.append(pointLevel_[pointi]);
    }


    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);
    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];
        cellMidPoint[celli] = 1234567890;    // mark need for splitting
    }


    boolList isDivisibleFace(mesh_.nFaces(), false);
    boolList isDivisibleEdge(mesh_.nEdges(), false);

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const label & patchID = mesh_.boundaryMesh().whichPatch(facei);

//        if (isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchID])
        if (mesh_.boundaryMesh()[patchID].type() == "empty")
        {
            isDivisibleFace[facei] = true;
            const labelList& fEdges = mesh_.faceEdges(facei);

            forAll(fEdges, i)
            {
                label edgeJ = fEdges[i];
                isDivisibleEdge[edgeJ] = true;
            }
        }
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels.size());

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                splitCells.insert(celli);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over cells to be refined or edges?
    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            const labelList& cEdges = mesh_.cellEdges(celli);
            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];
                const edge& e = mesh_.edges()[edgeI];
                if
                (
                 isDivisibleEdge[edgeI]
                 && pointLevel_[e[0]] <= cellLevel_[celli]
                 && pointLevel_[e[1]] <= cellLevel_[celli]
                )
                {
                    edgeMidPoint[edgeI] = 12345;    // mark need for splitting
                }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef4::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nFaces(); facei++)
    {
        faceAnchorLevel[facei] = getAnchorLevel(facei);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);


    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceAnchorLevel[facei] >= 0)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei];
            label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[facei]
             || newNeiLevel > faceAnchorLevel[facei]
            )
            {
                faceMidPoint[facei] = 1234567890;    // mark to be split (not a real nice way)
            }
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[facei] >= 0)
            {
                label own = mesh_.faceOwner()[facei];
                label ownLevel = cellLevel_[own];
                label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

                if
                (
                    newOwnLevel > faceAnchorLevel[facei]
                 || newNeiLevel[i] > faceAnchorLevel[facei]
                )
                {
                    faceMidPoint[facei] = 1234567890;    // mark to be split (not really nice way)
                }
            }
        }
    }


    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );



    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-GREAT, -GREAT, -GREAT)
        );

        forAll(bFaceMids, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceMidPoint[facei] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[facei];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0 && isDivisibleFace[facei])
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                const face& f = mesh_.faces()[facei];
                faceMidPoint[facei] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            facei < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[facei]
                          : bFaceMids[facei-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[facei]) = faceAnchorLevel[facei]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellLabels.size());

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                splitFaces.insert(facei);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split.
    // - faceMidPoint >= 0 : face needs to be split.
    // - edgeMidPoint >= 0 : edge needs to be split.
    // - isDivisibleFace true : face needs to be split in 4 else in 2.
    //   for the 2 split case, faceMidPoint is >= 0 but not initalized with a point.


    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());
    {
        labelList nAnchorPoints(mesh_.nCells(), 0);
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                cellAnchorPoints[celli].setSize(8);
            }
        }

       forAll(cellMidPoint, celli)
       {
            const cell& cFaces = mesh_.cells()[celli];
            forAll(cFaces, i)
            {
                label facei = cFaces[i];
                const face& f = mesh_.faces()[facei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if
                    (
                        cellMidPoint[celli] >= 0
                        && isDivisibleFace[facei]
                        && pointLevel_[pointi] <= cellLevel_[celli]
                    )
                    {
                        if (nAnchorPoints[celli] == 8)
                        {
                            dumpCell(celli);
                            FatalErrorIn
                                (
                                 "hexRef4::setRefinement(const labelList&"
                                 ", polyTopoChange&)"
                                ) << "cell " << celli
                                << " of level " << cellLevel_[celli]
                                << " uses more than 8 points of equal or"
                                << " lower level" << nl
                                << "Points so far:" << cellAnchorPoints[celli]
                                << abort(FatalError);
                        }
                        cellAnchorPoints[celli][nAnchorPoints[celli]++]
                            = pointi;
                    }
                }
            }
        }


        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                if (nAnchorPoints[celli] != 8)
                {
                    const labelList cPoints(cellPoints(celli));

                    FatalErrorIn
                    (
                        "hexRef4::setRefinement(const labelList&"
                        ", directTopoChange&)"
                    )   << "cell " << celli
                        << " of level " << cellLevel_[celli]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << IndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 3 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, celli)
    {
        const labelList& cAnchors = cellAnchorPoints[celli];

        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[celli];
            cAdded.setSize(4);

            // Original cell at 0
            cAdded[0] = celli;

            // Update cell level
            newCellLevel[celli] = cellLevel_[celli]+1;


            for (label i = 1; i < 4; i++)
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        celli,                              // master cell
                        mesh_.cellZones().whichZone(celli)  // zone for cell
                    )
                );

                newCellLevel(cAdded[i]) = cellLevel_[celli]+1;
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but only edges get split
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    affectedFace.set(cFaces[i]);
                }
            }
        }

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                affectedFace.set(facei);
            }
        }

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI);

                forAll(eFaces, i)
                {
                    affectedFace.set(eFaces[i]);
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, facei)
    {
        if (faceMidPoint[facei] >= 0 && affectedFace.get(facei))
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[facei];

            // Has original facei been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;
            label anchorLevel = faceAnchorLevel[facei];

            if(isDivisibleFace[facei])
            {
                face newFace(4);
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (pointLevel_[pointi] <= anchorLevel)
                    {
                        // point is anchor. Start collecting face.
                        DynamicList<label> faceVerts(4);
                        faceVerts.append(pointi);

                        // Walk forward to mid point.
                        // - if next is +2 midpoint is +1
                        // - if next is +1 it is midpoint
                        // - if next is +0 there has to be edgeMidPoint

                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        faceVerts.append(faceMidPoint[facei]);
                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        // Convert dynamiclist to face.
                        newFace.transfer(faceVerts);
                        // Get new owner/neighbour
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,          // Anchor point
                            own,
                            nei
                        );


                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];
                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }
                        if (!modifiedFace)
                        {
                            modifiedFace = true;
                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            }
            else
            {
                face newFace(2);
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    label nextpointi = f[f.fcIndex(fp)];
                    label edgeI = meshTools::findEdge (mesh_, pointi, nextpointi);
                    if (edgeMidPoint[edgeI] >=0)
                    {
                        DynamicList<label> faceVerts(4);
                        label pointJ = f[f.rcIndex(fp)];
                        faceVerts.append(pointi);
                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            f.rcIndex(fp),
                            faceVerts
                        );

                        faceVerts.append(pointJ);
                        newFace.transfer(faceVerts);
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,
                            own,
                            nei
                        );

                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];
                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }
                        if (!modifiedFace)
                        {
                            modifiedFace = true;
                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            }

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 2. faces that do not get split but use edges that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding edge splits to unsplit faces"
            << endl;
    }

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label facei = eFaces[i];

                if (faceMidPoint[facei] < 0 && affectedFace.get(facei))
                {
                    // Unsplit face. Add edge splits to face.

                    const face& f = mesh_.faces()[facei];
                    const labelList& fEdges = mesh_.faceEdges
                    (
                        facei,
                        fEdgesStorage
                    );

                    DynamicList<label> newFaceVerts(f.size());

                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] >= 0)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts);

                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    label anchorFp = findMinLevel(f);

                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        facei,
                        f[anchorFp],          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(facei))
                        {
                            label oldOwn = mesh_.faceOwner()[facei];
                            label oldNei = mesh_.faceNeighbour()[facei];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[facei];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[facei],
                                newFace
                            );
                        }
                    }

                    modFace(meshMod, facei, newFace, own, nei);

                    // Mark face as having been handled
                    affectedFace.unset(facei);
                }
            }
        }
    }


    // 3. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei))
        {
            const face& f = mesh_.faces()[facei];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                facei,
                f[anchorFp],          // Anchor point

                own,
                nei
            );

            modFace(meshMod, facei, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 4. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // This is the hard one. We have to find the splitting points between
    // the anchor points. But the edges between the anchor points might have
    // been split (into two,three or four edges).

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                celli,
                meshMod
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointi = labelMax;
        label maxPointi = labelMin;

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                minPointi = min(minPointi, cellMidPoint[celli]);
                maxPointi = max(maxPointi, cellMidPoint[celli]);
            }
        }
        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                minPointi = min(minPointi, faceMidPoint[facei]);
                maxPointi = max(maxPointi, faceMidPoint[facei]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointi = min(minPointi, edgeMidPoint[edgeI]);
                maxPointi = max(maxPointi, edgeMidPoint[edgeI]);
            }
        }

        if (minPointi != labelMax && minPointi != mesh_.nPoints())
        {
            FatalErrorIn("hexRef4::setRefinement(..)")
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointi:" << minPointi
                << " maxPointi:" << maxPointi
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexRef4::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, celli)
        {
            const labelList& addedCells = cellAddedCells[celli];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(celli, addedCells);
            }
        }
    }

    // Compact cellAddedCells.

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        refinedCells[i].transfer(cellAddedCells[celli]);
    }

    return refinedCells;
}


void Foam::hexRef4::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointi = pointsToStore[i];
        savedPointLevel_.insert(pointi, pointLevel_[pointi]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label celli = cellsToStore[i];
        savedCellLevel_.insert(celli, cellLevel_[celli]);
    }
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef4::updateMesh(const mapPolyMesh& map)
{
    Map<label> dummyMap(0);

    updateMesh(map, dummyMap, dummyMap, dummyMap);
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef4::updateMesh
(
    const mapPolyMesh& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef4::updateMesh :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef4::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef4 that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCellI)
            {
                label oldCellI = cellMap[newCellI];

                if (oldCellI == -1)
                {
                    newCellLevel[newCellI] = -1;
                }
                else
                {
                    newCellLevel[newCellI] = cellLevel_[oldCellI];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCellI = iter.key();
            label storedCellI = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCellI);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new cell "
                    << newCellI << " but cannot find old cell " << storedCellI
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCellI] = fnd();
        }

        //if (findIndex(cellLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef4::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "cellLevel_ contains illegal value -1 after mapping
        //        << " at cell " << findIndex(cellLevel_, -1) << endl
        //        << "This means that another program has inflated cells"
        //        << " (created cells out-of-nothing) and hence we don't know"
        //        << " their cell level. Continuing with illegal value."
        //        << abort(FatalError);
        //}
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef4 that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];

                if (oldPointi == -1)
                {
                    //FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    //    << "Problem : point " << newPointi
                    //    << " at " << mesh_.points()[newPointi]
                    //    << " does not originate from another point"
                    //    << " (i.e. is inflated)." << nl
                    //    << "Hence we cannot determine the new pointLevel"
                    //    << " for it." << abort(FatalError);
                    newPointLevel[newPointi] = -1;
                }
                else
                {
                    newPointLevel[newPointi] = pointLevel_[oldPointi];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointi = iter.key();
            label storedPointi = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointi);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new point "
                    << newPointi << " but cannot find old point "
                    << storedPointi
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointi] = fnd();
        }

        //if (findIndex(pointLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef4::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "pointLevel_ contains illegal value -1 after mapping"
        //        << " at point" << findIndex(pointLevel_, -1) << endl
        //        << "This means that another program has inflated points"
        //        << " (created points out-of-nothing) and hence we don't know"
        //        << " their point level. Continuing with illegal value."
        //        //<< abort(FatalError);
        //}
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.updateMesh(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.updateMesh(map);
}


// Gets called after mesh subsetting. Maps are from new back to old.
void Foam::hexRef4::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef4::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningIn
        (
            "hexRef4::subset(const labelList&, const labelList&"
            ", const labelList&)"
        )   << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCellI)
        {
            newCellLevel[newCellI] = cellLevel_[cellMap[newCellI]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef4::subset(..)")
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointi)
        {
            newPointLevel[newPointi] = pointLevel_[pointMap[newPointi]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef4::subset(..)")
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    //faceRemover_.subset(pointMap, faceMap, cellMap);
}


// Gets called after the mesh distribution
void Foam::hexRef4::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "hexRef4::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);
}


void Foam::hexRef4::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexRef4::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<> >
                    cellToFace(2*pp.size());

                label facei = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[facei];
                    label bFacei = facei-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFacei]), facei))
                    {
                        dumpCell(own);
                        FatalErrorIn("hexRef4::checkMesh()")
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << facei
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFacei]
                            << " has two faces connected to it:"
                            << facei << " and "
                            << cellToFace[labelPair(own, nei[bFacei])]
                            << abort(FatalError);
                    }

                    facei++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mag(mesh_.faceAreas()[i+mesh_.nInternalFaces()]);
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const scalar magArea = mag(mesh_.faceAreas()[facei]);

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[facei];
                label patchI = mesh_.boundaryMesh().whichPatch(facei);

                dumpCell(mesh_.faceOwner()[facei]);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[facei];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchI = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchI = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::checkMesh : Returning" << endl;
    }
}


void Foam::hexRef4::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexRef4::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn("hexRef4::checkRefinementLevels(const label)")
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << facei << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchI = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << facei
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronized
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointi)
        {
            if (pointLevel_[pointi] != syncPointLevel[pointi])
            {
                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointi << " coord:" << mesh_.points()[pointi]
                    << " has level " << pointLevel_[pointi]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointi]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointi)
        {
            const labelList& pCells = mesh_.pointCells(pointi);

            label& pLevel = maxPointLevel[pointi];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                label celli = pCells[i];

                if
                (
                    mag(cellLevel_[celli]-maxPointLevel[pointi])
                  > maxPointDiff
                )
                {
                    dumpCell(celli);

                    FatalErrorIn
                    (
                        "hexRef4::checkRefinementLevels(const label)"
                    )   << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << celli
                        << " cellLevel:" << cellLevel_[celli]
                        << " uses point:" << pointi
                        << " coord:" << mesh_.points()[pointi]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointi]
                        << abort(FatalError);
                }
            }
        }
    }
}



//
// Unrefinement
// ~~~~~~~~~~~~
//


////////////////

Foam::labelList Foam::hexRef4::getSplitEdges() const
{
    if (!history_.active())
    {
        FatalErrorIn("hexRef4::getSplitPoints()")<< "Only call if constructed with history capability"
            << abort(FatalError);
    }
    // Master cell
    // -1 undetermined
    // -2 certainly not split Edge
    // >= label of master cell
    labelList splitMaster(mesh_.nEdges(), -1);
    labelList splitMasterLevel(mesh_.nEdges(), 0);

    // Unmark all with not 4 cells (so internal edges)
    const labelListList& edgeCells = mesh_.edgeCells();

    forAll(edgeCells, edgeI)
    {
        const labelList& eCells = edgeCells[edgeI];

        if (eCells.size() != 4)
        {
            splitMaster[edgeI] = -2;
        }
    }
    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, celli)
    {
        const labelList& cEdges = mesh_.cellEdges()[celli];

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            label parentIndex = history_.parentIndex(celli);

            // Check same master.
            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                label masterCellI = splitMaster[edgeI];

                if (masterCellI == -1)
                {
                    // First time visit of edge. Store parent cell and
                    // level of the parent cell (with respect to celli). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[edgeI] = parentIndex;
                    splitMasterLevel[edgeI] = cellLevel_[celli] - 1;
                }
                else if (masterCellI == -2)
                {
                    // Already decided that point is not splitPoint
                }
                else if
                (
                    (masterCellI != parentIndex)
                 || (splitMasterLevel[edgeI] != cellLevel_[celli] - 1)
                )
                {
                    // Different masters so point is on two refinement
                    // patterns
                    splitMaster[edgeI] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                splitMaster[edgeI] = -2;
            }
        }
    }

    // Unmark boundary faces: maybe, it is probably a redundant check 
    // (maybe needed for parallel computations?)
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        const labelList& fEdges = mesh_.faceEdges()[facei];

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            splitMaster[edgeI] = -2;
        }
    }


    // Collect into labelList

    label nSplitEdges = 0;

    forAll(splitMaster, edgeI)
    {
        if (splitMaster[edgeI] >= 0)
        {
            nSplitEdges++;
        }
    }

    labelList splitEdges(nSplitEdges);
    nSplitEdges = 0;

    forAll(splitMaster, edgeI)
    {
        if (splitMaster[edgeI] >= 0)
        {
            splitEdges[nSplitEdges++] = edgeI;
        }
    }

    return splitEdges;
}


//void Foam::hexRef4::markIndex
//(
//    const label maxLevel,
//    const label level,
//    const label index,
//    const label markValue,
//    labelList& indexValues
//) const
//{
//    if (level < maxLevel && indexValues[index] == -1)
//    {
//        // Mark
//        indexValues[index] = markValue;
//
//        // Mark parent
//        const splitCell8& split = history_.splitCells()[index];
//
//        if (split.parent_ >= 0)
//        {
//            markIndex
//            (
//              maxLevel, level+1, split.parent_, markValue, indexValues);
//            )
//        }
//    }
//}
//
//
//// Get all cells which (down to level) originate from the same cell.
//// level=0 returns cell only, level=1 returns the 8 cells this cell
//// originates from, level=2 returns 64 cells etc.
//// If the cell does not originate from refinement returns just itself.
//void Foam::hexRef4::markCellClusters
//(
//    const label maxLevel,
//    labelList& cluster
//) const
//{
//    cluster.setSize(mesh_.nCells());
//    cluster = -1;
//
//    const DynamicList<splitCell8>& splitCells = history_.splitCells();
//
//    // Mark all splitCells down to level maxLevel with a cell originating from
//    // it.
//
//    labelList indexLevel(splitCells.size(), -1);
//
//    forAll(visibleCells, celli)
//    {
//        label index = visibleCells[celli];
//
//        if (index >= 0)
//        {
//            markIndex(maxLevel, 0, index, celli, indexLevel);
//        }
//    }
//
//    // Mark cells with splitCell
//}


Foam::labelList Foam::hexRef4::consistentUnrefinement
(
    const labelList& edgesToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef4::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorIn
        (
            "hexRef4::consistentUnrefinement(const labelList&, const bool"
        )   << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // Loop, modifying edgesToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect points to refine
    // maxSet = true: select points to refine

    // Maintain boolList for pointsToUnrefine and cellsToUnrefine
    PackedBoolList unrefineEdge(mesh_.nEdges());

    forAll(edgesToUnrefine, i)
    {
        label edgeI = edgesToUnrefine[i];

        unrefineEdge.set(edgeI);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const labelListList& edgeCells = mesh_.edgeCells();

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefineEdge, edgeI)
        {
            if (unrefineEdge.get(edgeI))
            {
                const labelList& pCells = edgeCells[edgeI];

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any point whose cell neighbour cannot be unrefined.
        forAll(unrefineEdge, edgeI)
        {
            if (unrefineEdge.get(edgeI))
            {
                const labelList& pCells = edgeCells[edgeI];

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefineEdge.unset(edgeI);
                        break;
                    }
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefineEdge, edgeI)
    {
        if (unrefineEdge.get(edgeI))
        {
            nSet++;
        }
    }

    labelList newEdgesToUnrefine(nSet);
    nSet = 0;

    forAll(unrefineEdge, edgeI)
    {
        if (unrefineEdge.get(edgeI))
        {
            newEdgesToUnrefine[nSet++] = edgeI;
        }
    }

    return newEdgesToUnrefine;
}


void Foam::hexRef4::setUnrefinement
(
    const labelList& splitEdgeLabels,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorIn
        (
            "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
        )   << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces(4*splitEdgeLabels.size());

        forAll(splitEdgeLabels, i)
        {
            const labelList& eFaces = mesh_.edgeFaces()[splitEdgeLabels[i]];

            forAll(eFaces, j)
            {
                splitFaces.insert(eFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorIn
            (
                "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "Ininitial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitEdgeLabels, i)
    {
        label edgeI = splitEdgeLabels[i];

        // Get original cell label

        const labelList& eCells = mesh_.edgeCells()[edgeI];

        // Check
        if (eCells.size() != 4)
        {
            FatalErrorIn
            (
                "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "splitEdge " << edgeI
                << " should have 4 cells using it. It has " << eCells
                << abort(FatalError);
        }


        // Check that the lowest numbered eCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCellI = min(eCells);

            forAll(eCells, j)
            {
                label celli = eCells[j];

                label region = cellRegion[celli];

                if (region == -1)
                {
                    FatalErrorIn("hexRef4::setUnrefinement(..)")
                        << "Ininitial set of split edges to unrefine does not"
                        << " seem to be consistent or not mid edges"
                        << " of refined cells" << nl
                        << "cell:" << celli << " on splitEdge " << edgeI
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCellI != cellRegionMaster[region])
                {
                    FatalErrorIn("hexRef4::setUnrefinement(..)")
                        << "cell:" << celli << " on splitEdge:" << edgeI
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the edgeCells:" << eCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the 4 cells that originated from merging around the split edge
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitEdgeLabels, i)
    {
        label edgeI = splitEdgeLabels[i];

        const labelList& eCells = mesh_.edgeCells()[edgeI];

        label masterCellI = min(eCells);

        forAll(eCells, j)
        {
            cellLevel_[eCells[j]]--;
        }

        history_.combineCells(masterCellI, eCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// Write refinement to polyMesh directory.
bool Foam::hexRef4::write() const
{
    bool writeOk = cellLevel_.write() && pointLevel_.write();

    if (history_.active())
    {
        writeOk = writeOk && history_.write();
    }

    return writeOk;
}


// ************************************************************************* //
