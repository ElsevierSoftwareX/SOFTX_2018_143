# README #

### Contributers ###

* Daniel Deising - Technical University Darmstadt
* Daniel Rettenmaier - Technical University Darmstadt
* T.G. Voskuilen ( https://github.com/tgvoskuilen/meshBalancing ) - Purdue University
* Timothée Pourpoint - Purdue University


### How do I get set up? ###
* `./Allwmake.sh`
* The `dynamic*Mesh` classes are a copy of the current **OF-4-dev** branch.
* If some solver needs dynamic linking of the libraries, make sure to source the new ones in `$FOAM_USER_LIBBIN`
* Use the `dynamicMeshDict` as provided in the folder `dynamicRefineBalancedFvMesh`.
* Use the test cases as a basis for new ones.

### Who do I talk to? ###

* Daniel Rettenmaier: rettenmaier@gsc.tu-darmstadt.de
* Daniel Deising: deising@mma.tu-darmstadt.de

***
### Utilities ###
###### meshUpdater
Updates the mesh calling `mesh.update()`. Depending on the selected mesh type in `dynamicMeshDict`, an adaptive mesh refinement and a load balancing step is performed.
Thereby all volume and surface fields are mapped onto the new mesh or processor distribution. 

An overwrite option is provided using `-overwrite`

###### initSurfaceFields
Initializes a `surfaceScalarField mySurfaceScalarField;` and a `surfaceVectorField mySurfaceVectorField;` which are initialized with the distance to the domain origin.


### Test Cases ###
###### unitTest_surfaceFieldMapping
A domain with two cells is initialized using the utility `initSurfaceFields`. Calling `meshUpdater` the domain is refined in one step into 16 cells. This creates 12 internal patches which are of main interest here. Our library setup tackles two serious issues which have been identified using this case: 
*    Wrong addressing of internal patches: This problem arises in the function call [`hexRef8::addInternalFace()`](https://bitbucket.org/drettenmaier/amrandloadbalancing/commits/1cef6f9569fe7525b5b769c43c910ec50dbd6785?at=master#chg-src/dynamicMesh/polyTopoChange/polyTopoChange/hexRef8/hexRef8.C), where a new internal face is either "expanded" from a master face or without a master face. As we understand it, the first option makes no sense and furthermore changes the face index to meshFacei, which in the following corresponds to some other face which is not related with the new one. The newly created internal face always need a face index of -1. 

*    No mapping of non flux surfaceFields:

###### unitTest_LoadBalancingSurfaceFieldFlip

###### damBreakAMR_LB

### Code ###