
# README #

[TOC]

* * *

### Value proposition
Dynamic cell refinement reduces computational costs for problems that changes over time by adapting the cell size automatically to accuracy requirements. For cases that are calculated in parallel the load between processors might become imbalanced resulting in idle processor time. 

The dynamic refinement algorithm of OpenFOAM is enhanced by  

* adding very flexible refinement based on multiple criteria 
* and additional to the 3D a 2D and 2.5D refinement
* fixing severe bugs regarding 
    - face addressing, 
    - mapping of cell and boundary fields between refinement states
    - and flux flipping

Load balancing is made possible by a sophisticated redistribution algorithm already implemented in OpenFOAM. The load balancing code of T.G. Voskuilen is used and stabilized fixing bugs regarding

* flux flipping of non flux fields
* and boundary mapping 

Load balancing is used in the test cases exemplarily for, but not limited to, the interDyMFoam solver. 

Two unit cases illustrate the main bug fixes and several tutorial cases show possible use cases in the area of multiphase flows of dynamic mesh refinement in combination with load balancing.

* * *

### Contributors ###

* Daniel Rettenmaier - Technical University Darmstadt
* Daniel Deising - Technical University Darmstadt
* Yun Ouedraogo - Technical University Darmstadt
* Holger Marschall - Technical University Darmstadt
* Other (indirect) contributors
    - T.G. Voskuilen ([meshBalancing](https://github.com/tgvoskuilen/meshBalancing)) - Purdue University (Dynamic load balancing in OpenFOAM-2.3.x)
    - Timothee Pourpoint - Purdue Universit (Dynamic load balancing in OpenFOAM-2.3.x)  
    - Ahmad Baniabedalruhman - Michigan Technological University ([2D refinement](http://faculty.yu.edu.jo/ahmad_a/Lists/Other%20Academic%20Activities/AllItems.aspx))
    - Stefan Batzdorf - Technical University Darmstadt ([`decomposeParLevel`](https://bitbucket.org/drettenmaier/amr_lb_publication/src/master/OpenFOAM/applications/utilities/reconstructParLevel/reconstructParLevel.C))
    - Andrea Montorfano - Politecnico de Milano  (Discussions and hint for dontFlipPatch)

* * *

### How do I get set up? ###
* Checkout the branch which fits your OpenFOAM installation
    - v-dev (commit OpenFOAM-dev: 9dcdf23a6b94a8d792d94664ccfd0d7948a5c905 )
    - v5.x (commit OpenFOAM-dev: d06c3b390ac18dc2435ae87330038265a69c1c56 )
    - v4.x (commit OpenFOAM-dev: d214c8dfd5ba56dd442bae186fd4fb50dd35c338 )
* Copy a patch into your OpenFOAM installation path
  ```cp -r OpenFOAM/dontFlipSurfaceVectorFields.patch $WM_PROJECT_DIR```
* Have a short look inside the `*.patch` file. It is a minor change to your source code!
* Apply `git apply dontFlipSurfaceVectorFields.patch` in OpenFOAM core installation. 
* Compile OpenFOAM core again
* Source project specific variables ```source OpenFOAM/etc/bashrc```
* Compile libraries and applications of this repositories `./Allwmake.sh`
* If some solver needs dynamic linking of the libraries, make sure to source the new ones in `$FOAM_USER_LIBBIN`
* Take note that the corresponding libraries are called `dynamic*Mesh.so`
* Use the `dynamicMeshDict` as provided in here `src/dynamicRefineBalancedFvMesh/dynamicMeshDict`.

### Who do I talk to? ###

* Daniel Rettenmaier: rettenmaier@gsc.tu-darmstadt.de
* Daniel Deising: deising@mma.tu-darmstadt.de
* Yun Ouedraogo: ouedraogo@temf.tu-darmstadt.de 
* Holger Marschall: marschall@mma.tu-darmstadt.de 

***

### Utilities ###
###### meshUpdater
Updates the mesh calling `mesh.update()`. Depending on the selected mesh type in `dynamicMeshDict`, an adaptive mesh refinement and a load balancing step is performed.
Thereby all volume and surface fields are mapped onto the new mesh or processor distribution. 

An overwrite option is available `-overwrite`

###### meshUpdaterOrig
Does the same as meshUpdater, it is however linked with the old `dynamicFvMesh` to compare the difference in the unit test cases.

###### initSurfaceFields
A utility used for unit test cases which initializes a `surfaceScalarField mySurfaceScalarField;` and a `surfaceVectorField mySurfaceVectorField;` which are initialized with the distance to the domain origin.

###### decomposeParLevel
A development from Christian Kunkelmann and Stefan Batzdorf (2015) to ensure a correct decomposition of the refinement history and the refinement levels. 
This utility is only necessary in 2D and 2.5D cases. In 3D the refinement constraints in the `decomposeParDict` or `balanceParDict` will do the trick.

###### reconstructParLevel
A development from Christian Kunkelmann and Stefan Batzdorf (2015) to ensure a correct reconstruction of the refinement history and the refinement levels. This utility is only necessary in 2D and 2.5D cases. In 3D the refinement constraints in the `decomposeParDict` or `balanceParDict` will do the trick.

* * *

### Test Cases ###
##### surfaceFieldMapping - unitTest
A domain with two cells is initialized using the utility `initSurfaceFields`. Calling `meshUpdater` the domain is refined in one step into 16 cells. This creates 12 internal patches which are of main interest here. Our library setup tackles two serious issues which have been identified using this case: 

*    **Wrong addressing of internal patches**: This problem arises in the function call [`hexRef::addInternalFace()`](https://bitbucket.org/drettenmaier/amr_lb_publication/src/v-dev/OpenFOAM/src/dynamicMesh/polyTopoChange/polyTopoChange/hexRef/hexRef.C), where a new internal face is either "expanded" from a master face or without a master face. As we understand it, the first option makes no sense and furthermore changes the face index to meshFacei, which in the following corresponds to some other face which is not related with the new one. The newly created internal face always need a face index of -1. So we simply removed the option of creating an internal face with a corresponding master-face index.

*    **No mapping of non flux surfaceFields**: Currently the newly created internal faces have no build in mapping mechanism of surfaceFields at the refinement step. Only fluxes, which are listed in `dynamicMeshDict`, are recalculated after the mesh refinement step, which only works for surfaceScalarFields. So we introduced a surfaceFields mapping for fields which are listed in a new list called mapSurfaceFields in [`dynamicMeshDict`](https://bitbucket.org/drettenmaier/amr_lb_publication/src/a32f2ac2e862694962dcd6751467e9925193e715/OpenFOAM/src/dynamicFvMesh/dynamicRefineBalancedFvMesh/dynamicMeshDict?at=master&fileviewer=file-view-default). The fields of internal faces are averaged using the three adjacent unrefined face values.

##### surfaceFieldMapping2D - unitTest
Analog to the 3D case. However, the front and back boundary patches are empty and therefore the surfaceFields have a default value on the empty patches.

##### LoadBalancingSurfaceFieldFlip - unitTest
A domain with three cells is decomposed manually. processor0 with one and processor1 with two cells. The cell at processor0 gets refined which leads to an imbalance higher than the minimal threshold for load balancing. **During the rebalancing step all surfaceFields get mapped using the flipping mechanism** as it is only correct for fluxes. Signs of the surfaceFields are falsely changed at processor patches if they don't represent a flux.
Code changes are were necessary in [`fvMeshDistribute/fvMeshDistributeTemplates.C`](https://bitbucket.org/drettenmaier/amr_lb_publication/commits/1cef6f9569fe7525b5b769c43c910ec50dbd6785?at=master#chg-src/dynamicMesh/fvMeshDistribute/fvMeshDistributeTemplates.C), [`fvMeshSubset/fvMeshSubsetInterpolate.C`](https://bitbucket.org/drettenmaier/amr_lb_publication/commits/1cef6f9569fe7525b5b769c43c910ec50dbd6785?at=master#chg-src/dynamicMesh/fvMeshSubset/fvMeshSubsetInterpolate.C) and [`src/finiteVolume/interpolation/mapping/fvFieldMappers/MapFvSurfaceField.H`](https://bitbucket.org/drettenmaier/amr_lb_publication/commits/1cef6f9569fe7525b5b769c43c910ec50dbd6785?at=master#chg-dontFlipSurfaceVectorFields.patch)
In some cases we still observe a flipping of `surfaceVectorFields` using the `decomposePar` method. The work-flow of initializing a surfaceVectorField and then decompose the domain is therefore not recommended.

##### damBreakAMR_LB
The damBreak case is set up with an adaptive mesh refinement and load balancing. `interDyMFoam` recalculates the flux `phi` after a mesh.update() call using the surfaceVectorField `Uf`, where in earlier versions of OF `fvc::interpolate(U)` was used. Without the appropriate mapping of new internal faces and the careful handling of the face flipping operator AMR+LB simulations only run through with some luck.

##### capillaryRisePlate2D and in a Tube 2.5D
CapillaryRise between two plates

Since the meniscus resembles the only point of interest a local quadtree refinement reduces computational effort significantly by keeping a high resolution. The meniscus rises due to the contact angle and therefore load balancing is useful to avoid bottlenecks.

A similar test case is the capillary rise in a tube where the axisymmetric 2.5D refinement is applied.

When changing the mesh OpenFOAM provides the functionality to map a value between those changes. A the index of a later refined cell is mapped to its child cells and vice versa. In the current OpenFOAM versions however no default mapping is implemented for boundaries that hold a gradient, such as `alphaContactAngleFvPatchField.H` and its base-class `fixedGradientFvPatchField.H` or `mixedFvPatchField.H`. The gradient representation is reallocated in size but not properly initialized. A call for `evaluate()` on the boundary at the load balancing might therefore lead to an error. 

This issue is fixed exemplarily for the [alphaContactAngle](https://bitbucket.org/drettenmaier/amr_lb_publication/src/a32f2ac2e862694962dcd6751467e9925193e715/OpenFOAM/src/twoPhaseProperties/alphaContactAngle/alphaContactAngleFvPatchScalarField.H?at=master&fileviewer=file-view-default) by specifying the mapper. The gradient on newly added faces is set to zero. To be save make sure that cells containing the three-phase contact line are not refined. Note that when using e.g. `fixedFluxPressure` a similar specification is necessary. 

* * *

### Scripts ###

The folder `OpenFOAM/scripts/` should be sourced with `source OpenFOAM/etc/bashrc` to make all scripts available wherever your test case is located. 

##### reconstructParMeshAll.sh

If the mesh changes a standard `reconstructPar` will not be sufficient to rearrange the mesh. Therefore `reconstructParMesh` needs to be called first. And for 2D and 2.5D cases `reconstructParLevel` assures that the refinement history gets its proper treatment. 
All those steps are combined in `OpenFOAM/scripts/reconstructParMeshAll.sh`. It automatically checks for all time steps that need to be reconstructed and will do all necessary steps. 
Options like `-time '0:1.2'`,  `-latestTime`  or `-fields '(alpha.water)'` are often helpful and will be piped through the bash script to work accordingly.

##### cleanCase.sh

Will remove all time steps, even your `0` folder, all processor folders, `log.*` files and the mesh in `constant/polyMesh`. This script is used to cleanup the tutorial cases. 

##### touchCaseFoam.sh

Will create a `case.foam` file in your case and in each `processor*` folder so you might load sub-domains separately into Paraview. 


