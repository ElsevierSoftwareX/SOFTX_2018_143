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
### utilities ###
###### meshUpdater
Updates the mesh calling `mesh.update()`. Depending on the selected mesh type in `dynamicMeshDict`, an adaptive mesh refinement and a load balancing step is performed.
Thereby all volume and surface fields are mapped onto the new mesh or processor distribution. 

An overwrite option is provided using `-overwrite`

###### initSurfaceFields
Initializes a `surfaceScalarField mySurfaceScalarField;` and a `surfaceVectorField mySurfaceVectorField;` which are initialized with the distance to the domain origin.


### Test Cases ###

### Code ###