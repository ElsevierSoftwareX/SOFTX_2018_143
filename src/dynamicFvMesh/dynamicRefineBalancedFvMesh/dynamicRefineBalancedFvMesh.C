/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "dynamicRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineBalancedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

Foam::List<Foam::scalar> Foam::dynamicRefineBalancedFvMesh::readRefinementPoints()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );
    
    List<scalar> refData(4, scalar(0));
    
    refData[0] = readScalar(refineDict.lookup("unrefineLevel"));
    refData[1] = readScalar(refineDict.lookup("lowerRefineLevel"));
    refData[2] = readScalar(refineDict.lookup("refineInterval"));
    refData[3] = readScalar(refineDict.lookup("maxRefinement"));
        
    return refData;
}

void Foam::dynamicRefineBalancedFvMesh::updateRefinementField()
{
    Info<< "Calculating internal refinement field" << endl;
    
    volScalarField& intRefFld = *internalRefinementFieldPtr_;
    volScalarField& targetFld = *targetLevelPtr_;
    volScalarField& currFld = *isLevelPtr_;
    
    // Set the internal refinement field to zero to start with
    intRefFld = dimensionedScalar("zero",dimless,0.0);
    
    // Get the cell level field from dynamicRefineFvMesh
    const labelList& cellLevel = meshCutter().cellLevel();
    
    // Read the points at which refinement and unrefinement occur from the
    // dynamicMeshDict entries
    List<scalar> refinePoints = readRefinementPoints();
    
    // init list for refine/unrefine field (-1=unrefine, 1=refine, 0=do nothing)
    labelList markRefineCells (this->nCells(), 0);
    
    // init list for target refinement level per cell
    labelList targetLevel (this->nCells(), 0);
    
    // First fields
    List<word> fieldNames = fields_.toc();
    Field<scalar> refFld(nCells(),0.0);
    
    forAll(fieldNames, i)
    {
        word fldName = fieldNames[i];

        scalar minValue = readScalar(fields_[fldName].lookup("minValue"));
        scalar maxValue = readScalar(fields_[fldName].lookup("maxValue"));
        scalar refineLevel = readScalar(fields_[fldName].lookup("refineLevel"));
        
        const volScalarField& fld = this->lookupObject<volScalarField>(fldName);
        
        // Limit the value of refFld based on its max level
        forAll(fld, cellI)
        {
            if ((fld[cellI] >= minValue) && (fld[cellI] <= maxValue))
            {
                // increase targetLevel up to refineLevel 
                // BUT: do not decrease if cell already marked for higher refinement level by previous criterion
                targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
            }
        }
        
        //- AddLayer Keyword
        scalar nAddLayers(0);
        if (fields_[fldName].found("nAddLayers"))
        {
            nAddLayers = readScalar(fields_[fldName].lookup("nAddLayers"));
        }

/*
        Info << "Foam::dynamicRefineBalancedFvMesh::updateRefinementField FieldName: " <<fldName <<endl;
        Info << "Foam::dynamicRefineBalancedFvMesh::updateRefinementField minValue: " <<minValue <<endl;
        Info << "Foam::dynamicRefineBalancedFvMesh::updateRefinementField maxValue: " <<maxValue <<endl;
        Info << "Foam::dynamicRefineBalancedFvMesh::updateRefinementField refineLevel: " <<refineLevel <<endl;
        Info << "Foam::dynamicRefineBalancedFvMesh::updateRefinementField addLayers: " <<nAddLayers <<endl;
*/

        if (nAddLayers > 0)
        {
            //- create a volField of the labelList targetLevel
            volScalarField tLevel(intRefFld * 0.0);  
            forAll(targetLevel, cellI)
            {
                tLevel[cellI] = targetLevel[cellI];         
            }

            // #nAddLayers cells per additional level
            for(label j=0; j < nAddLayers; j++)
            {
                //- select the area with targetLevel==refineLevel
                volScalarField finest = pos(tLevel - refineLevel + SMALL);
                
                //- add +1 to targetLevel on the enlarged stencil 
                tLevel += pos( fvc::average(fvc::interpolate(finest) - SMALL)) - finest;
            }

            //- copy the new results onto the target labelList
            forAll(targetLevel, cellI)
            {
                targetLevel[cellI] = tLevel[cellI];         
            }
        }
    }
    
    // Then gradients
    List<word> gradFieldNames = gradFields_.toc();
    
    Field<scalar> cubeRtV = Foam::pow(this->V(),1.0/3.0);
    
    forAll(gradFieldNames, i)
    {
        word fldName = gradFieldNames[i];
        scalar minValue = gradFields_[fldName][0];
        scalar maxValue = gradFields_[fldName][1];
        label refineLevel = static_cast<label>(gradFields_[fldName][2]);
        
        const volScalarField& fld = this->lookupObject<volScalarField>(fldName);

        refFld = mag(fvc::grad(fld)) * cubeRtV;
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if ((refFld[cellI] >= minValue) && (refFld[cellI] <= maxValue))
            {
                targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
            }
        }
    }
    
    // Then curls
    List<word> curlFieldNames = curlFields_.toc();
    
    
    forAll(curlFieldNames, i)
    {
        word fldName = curlFieldNames[i];
        scalar minValue = curlFields_[fldName][0];
        scalar maxValue = curlFields_[fldName][1];
        label refineLevel = static_cast<label>(curlFields_[fldName][2]);
        
        const volVectorField& fld = this->lookupObject<volVectorField>(fldName);
        
        refFld = mag(fvc::curl(fld));
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if ((refFld[cellI] >= minValue) && (refFld[cellI] <= maxValue))
            {
                targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
            }
        }
    }
    
    // At the interface, assumed to be always the maximum refinement level
    List<word> interfaceRefineField = interface_.toc();
    forAll(interfaceRefineField, i)
    {
        word fldName = interfaceRefineField[i];
        // read region of maximum refinement levels inside and outside of interface indicator field 
        // (does not need to be alpha, can also be a concentration field)
        
        scalar innerRefLayers = readScalar(interface_[fldName].lookup("innerRefLayers"));
        scalar outerRefLayers = readScalar(interface_[fldName].lookup("outerRefLayers"));
        
        scalar nAddLayers(0);
        if (interface_[fldName].found("nAddLayers"))
        {
            nAddLayers = readScalar(interface_[fldName].lookup("nAddLayers"));
        }
        
        label refineLevel = refinePoints[3];
        
        if (interface_[fldName].found("maxRefineLevel"))
        {
            refineLevel = readScalar(interface_[fldName].lookup("maxRefineLevel"));

            // to avoid user input mistakes, limit the value with the maximum allowed
            refineLevel = min(refinePoints[3], refineLevel);
        }
        
        const volScalarField& fld = this->lookupObject<volScalarField>(fldName);
        
        volScalarField isInterface(intRefFld * 0.0);
        isInterface = dimensionedScalar("isInterface",dimless,0.0);
        
        surfaceScalarField deltaAlpha = mag(fvc::snGrad(fld) / this->deltaCoeffs());
        
        const unallocLabelList& owner = this->owner();
        const unallocLabelList& neighbour = this->neighbour();
        
        forAll(deltaAlpha, faceI)
        {
        //TODO: add a min max value of the specified field to be marked for refinement
        // NOW: hard-coded method snGradAlpha: checks mag(alpha_N - alpha_P) > 0.1 ? 1:0
            if (deltaAlpha[faceI] > 0.1) // currently a fixed prescribed value; should be read from díctionary
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];
                
                // set isInterface field to one
                isInterface[own] = 1.0;
                isInterface[nei] = 1.0;
            }
        }
        
        // assumed fld=0.5*(fldMax+fldMin) defines the interface
        dimensionedScalar fldInterfaceValue(0.5*(max(fld)+min(fld)));
        
        //-DD: old implementation based on face interpolation
        //     which results in slower transport in diagonal direction
        // add inner refinement layers
        //for(label i=0; i < innerRefLayers; i++)
        //{
        //    isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos(fld - fldInterfaceValue)); 
        //    isInterface = neg(- isInterface);
        //}
        //
        // add outer refinement layers
        //for(label i=0; i < outerRefLayers; i++)
        //{
        //    isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos(fldInterfaceValue - fld));
        //    isInterface = neg(- isInterface);
        //}
        //
        //forAll(isInterface, cellI)
        //{
        //    if (isInterface[cellI] > 0.5)
        //    {
        //        targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
        //    }
        //}

        //-DD: new version using volPointInterpolation (direction independent buffer layer)
        const volPointInterpolation& pInterp = volPointInterpolation::New(*this);
        const fvMesh& mesh = fld.mesh();
     
        // add inner refinement layers
        for(label i=0; i < innerRefLayers; i++)
        {
            volScalarField markInner(isInterface*pos(fld - fldInterfaceValue));
            pointScalarField markLayerP(pInterp.interpolate(markInner));

            forAll(mesh.C(), cellI)
            {
                scalar sum = 0.;
                label nPoints = 0;

                forAll(mesh.cellPoints()[cellI], pointI)
                {
                    sum += markLayerP[mesh.cellPoints()[cellI][pointI]];
                    nPoints++;
                }
                if (nPoints > 0)
                {
                    sum /= nPoints;
                }
                isInterface[cellI] += sum;
            }
        }
        isInterface = pos(isInterface - SMALL);

        // add outer refinement layers
        for(label i=0; i < outerRefLayers; i++)
        {
            volScalarField markOuter(isInterface*pos(fldInterfaceValue - fld));
            pointScalarField markLayerP(pInterp.interpolate(markOuter));

            forAll(mesh.C(), cellI)
            {
                scalar sum = 0.;
                label nPoints = 0;

                forAll(mesh.cellPoints()[cellI], pointI)
                {
                    sum += markLayerP[mesh.cellPoints()[cellI][pointI]];
                    nPoints++;
                }
                if (nPoints > 0)
                {
                    sum /= nPoints;
                }
                isInterface[cellI] += sum;
            }
        }
        isInterface = pos(isInterface - SMALL);

        forAll(isInterface, cellI)
        {
            if (isInterface[cellI] > 0.5)
            {
                targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
            }
        }
        
        //-DD: old implementation based on face interpolation
        //     which results in slower transport in diagonal direction
        // expand additional layers if specified:
        if (nAddLayers > 0)
        {
           // loop over additional layers
           for(label i=1; i < refineLevel; i++)
           {
               // #nAddLayers cells per additional level
               for(label j=0; j < nAddLayers; j++)
               {
                   isInterface += neg(- fvc::average(fvc::interpolate(isInterface)));
                   isInterface = neg(- isInterface);
               }
               
               forAll(isInterface, cellI)
               {
                   if (isInterface[cellI] == 1.0)
                   {
                       targetLevel[cellI] = max(targetLevel[cellI], (refineLevel - i));
                   }
               }
           }
        }
    }
    
    // The set refinement physical regions (force the mesh to stay refined
    // near key features)
    forAll(refinedRegions_, regionI)
    {
        const entry& region = refinedRegions_[regionI];
        
        autoPtr<topoSetSource> source = 
            topoSetSource::New(region.keyword(), *this, region.dict());
        
        cellSet selectedCellSet
        (
            *this,
            "cellSet",
            nCells()/10+1 //estimate
        );
        
        source->applyToSet
        (
            topoSetSource::NEW,
            selectedCellSet
        );
        
        const labelList cells = selectedCellSet.toc();
        
        label minLevel = readLabel(region.dict().lookup("minLevel"));
        
        forAll(cells, i)
        {
            const label& cellI = cells[i];
            
            targetLevel[cellI] = max(targetLevel[cellI], minLevel);
        }
    }
    
    //-DD: add buffer layer based on targetLevel field to prevent 2-to-1 refinement
    
    volScalarField blockedLevel = targetFld * 0.;
    
    for (label currLayer=refinePoints[3]; currLayer>=1; currLayer--)
    {
        forAll (targetLevel, cellI)
        {
            if (targetLevel[cellI] >= currLayer)
            {
                blockedLevel[cellI] = targetLevel[cellI];
            }
        }
        
        for (label i=0; i<nBufferLayers_; i++)
        {
            blockedLevel = max(blockedLevel, pos(fvc::average(fvc::interpolate(blockedLevel)) - SMALL)*(currLayer-1));
        }

        labelList blockLev(blockedLevel.internalField().size(), 0);
        forAll (blockLev, cellI)
        {
            blockLev[cellI] = blockedLevel[cellI];
        }
        targetLevel = max(targetLevel, blockLev);
    }
    
    // mark cells to be refined/unrefined based on above criteria:
    // simply check, if targetLevel lower or higher cellLevel
    forAll (intRefFld.internalField(), cellI)
    {
        intRefFld.primitiveFieldRef()[cellI] = targetLevel[cellI] - cellLevel[cellI];
        targetFld.primitiveFieldRef()[cellI] = targetLevel[cellI];
        currFld.primitiveFieldRef()[cellI] = cellLevel[cellI];
    }

    intRefFld.correctBoundaryConditions();
    
    Info<<"Min,max refinement field = " << Foam::min(intRefFld).value() << ", "
        << Foam::max(intRefFld).value() << endl;
}

void Foam::dynamicRefineBalancedFvMesh::readRefinementDict()
{
    dictionary dynamicMeshDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    nBufferLayers_ = readLabel(dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs").lookup("nBufferLayers"));
    
    if( dynamicMeshDict.isDict("refinementControls") )
    {
        dictionary refineControlDict = 
            dynamicMeshDict.subDict("refinementControls");
        
        enableRefinementControl_ = 
            Switch(refineControlDict.lookup("enableRefinementControl"));
        
        if( enableRefinementControl_ )
        {
            // Overwrite field name entry in dynamicRefineFvMeshCoeffs?
            // For now you just have to be smart and enter 
            // 'internalRefinementField' for the field name manually
            
            // Read HashTable of field-refinement scalars
            if( refineControlDict.found("fields") )
            {
                fields_ = HashTable< dictionary >
                (
                    refineControlDict.lookup("fields")
                );
            }
            
            // Read HashTable of gradient-refinement scalars
            if( refineControlDict.found("gradients") )
            {
                gradFields_ = HashTable< List<scalar> >
                (
                    refineControlDict.lookup("gradients")
                );
            }
            
            // Read HashTable of curl-refinement vectors
            if( refineControlDict.found("curls") )
            {
                curlFields_ = HashTable< List<scalar> >
                (
                    refineControlDict.lookup("curls")
                );
            }
            
            // Read interface refinement data
            if( refineControlDict.found("interface") )
            {
                interface_ = HashTable< dictionary >
                (
                    refineControlDict.lookup("interface")
                );
            }
            
            // Read refinement regions
            if( refineControlDict.found("regions") )
            {
                refinedRegions_ = PtrList<entry>
                (
                    refineControlDict.lookup("regions")
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
            
Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io),
    internalRefinementFieldPtr_(NULL),
    targetLevelPtr_(NULL),
    isLevelPtr_(NULL),
    fields_(),
    gradFields_(),
    curlFields_(),
    interface_(),
    refinedRegions_(),
    enableRefinementControl_(false),
    nBufferLayers_(0),
    rebalance_(false)
{
    readRefinementDict();
    
    if( enableRefinementControl_ )
    {
        internalRefinementFieldPtr_ = new volScalarField
        (
            IOobject
            (
                "internalRefinementField",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );

        targetLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "targetLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );

        isLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "isLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::update()
{
    //Part 0 - Update internally calculated refinement field
    readRefinementDict();
    
    List<scalar> refinePoints = readRefinementPoints();
    label refineInterval = refinePoints[2];
    
    if( enableRefinementControl_ && time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0 )
    {
        updateRefinementField();
    }

    //Part 1 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();

    if( Pstream::parRun() && hasChanged)
    {
        //Correct values on all coupled patches
        correctBoundaries<scalar>();
        correctBoundaries<vector>();
        correctBoundaries<sphericalTensor>();
        correctBoundaries<symmTensor>();
        correctBoundaries<tensor>();      
    }

    dictionary decomposeParDict
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                time().system(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    
    rebalance_ = false; 

    // Part 2 - Load Balancing
    {    
        dictionary refineDict
        (
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    time().constant(),
                    *this,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).subDict("dynamicRefineFvMeshCoeffs")
        );
    
        Switch enableBalancing = refineDict.lookup("enableBalancing");
    
        if ( Pstream::parRun() && hasChanged )
        {
            const scalar allowableImbalance = 
                readScalar(refineDict.lookup("allowableImbalance"));
            
            //First determine current level of imbalance - do this for all
            // parallel runs with a changing mesh, even if balancing is disabled
            label nGlobalCells = globalData().nTotalCells();
            scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
            scalar localImbalance = mag(scalar(nCells()) - idealNCells);
            Foam::reduce(localImbalance, maxOp<scalar>());
            scalar maxImbalance = localImbalance/idealNCells;
        
            Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;
        
            //If imbalanced, construct weighted coarse graph (level 0) with node
            // weights equal to their number of subcells. This partitioning works
            // as long as the number of level 0 cells is several times greater than
            // the number of processors.
            if( maxImbalance > allowableImbalance && enableBalancing)
            {
                Info << "\n**Solver hold for redistribution at time = "  << time().timeName() << " s" << endl;
                
                rebalance_ = true;   

                const labelIOList& cellLevel = meshCutter().cellLevel();
                Map<label> coarseIDmap(100);
                labelList uniqueIndex(nCells(),0);
            
                label nCoarse = 0;

                forAll(cells(), cellI)
                {
                    if( cellLevel[cellI] > 0 )
                    {
                        uniqueIndex[cellI] = nCells() + topParentID
                        (
                            meshCutter().history().parentIndex(cellI)
                        );
                    }
                    else
                    {
                        uniqueIndex[cellI] = cellI;
                    }
                
                    if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                    {
                        ++nCoarse;
                    }
                }
            
                // Convert to local sequential indexing and calculate coarse
                // points and weights
                labelList localIndex(nCells(),0);
                pointField coarsePoints(nCoarse,vector::zero);
                scalarField coarseWeights(nCoarse,0.0);
            
                forAll(uniqueIndex, cellI)
                {
                    localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];
                
                    // If 2D refinement (quadtree) is ever implemented, this '3'
                    // should be set in general as the number of refinement
                    // dimensions.
                    label w = (1 << (3*cellLevel[cellI]));
                
                    coarseWeights[localIndex[cellI]] += 1.0;
                    coarsePoints[localIndex[cellI]] += C()[cellI]/w;
                }
            
                // Set up decomposer - a separate dictionary is used here so
                // you can use a simple partitioning for decomposePar and
                // ptscotch for the rebalancing (or any chosen algorithms)
                autoPtr<decompositionMethod> decomposer
                (
                    decompositionMethod::New
                    (
                        IOdictionary
                        (
                            IOobject
                            (
                                "balanceParDict",
                                time().system(),
                                *this,
                                IOobject::MUST_READ_IF_MODIFIED,
                                IOobject::NO_WRITE
                            )
                        )
                    )
                );
            
                labelList finalDecomp = decomposer().decompose
                (
                    *this, 
                    localIndex,
                    coarsePoints,
                    coarseWeights
                );

                scalar tolDim = globalMeshData::matchTol_ * bounds().mag();
            
                Info<< "Distributing the mesh ..." << endl;
                fvMeshDistribute distributor(*this, tolDim);
            
                Info<< "Mapping the fields ..." << endl;
                autoPtr<mapDistributePolyMesh> map =
                    distributor.distribute(finalDecomp);
                  
                Info<< "Distribute the map ..." << endl;
                meshCutter_.distribute(map);


                Info << "Successfully distributed mesh" << endl;

                scalarList procLoadNew (Pstream::nProcs(), 0.0);
                procLoadNew[Pstream::myProcNo()] = this->nCells();

                reduce(procLoadNew, sumOp<List<scalar> >());

                scalar overallLoadNew = sum(procLoadNew);
                scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

                Info << "New distribution: " << procLoadNew << endl;
                Info << "Max deviation: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;
            }
        }
    }

    if( Pstream::parRun() && rebalance_)
    {
        //Correct values on all coupled patches
        correctBoundaries<scalar>();
        correctBoundaries<vector>();
        correctBoundaries<sphericalTensor>();
        correctBoundaries<symmTensor>();
        correctBoundaries<tensor>();      
    }

    return hasChanged;
}


// ************************************************************************* //
