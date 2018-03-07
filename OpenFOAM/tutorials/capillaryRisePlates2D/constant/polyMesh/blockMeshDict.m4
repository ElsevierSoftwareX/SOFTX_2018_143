/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(int_, [esyscmd(perl -e 'print ($1)')]) dnl>

define(hex2D, hex (b$1 b$2 b$3 b$4 f$1 f$2 f$3 f$4))
define(quad2D, (b$1 b$2 f$2 f$1))
define(frontQuad, (f$1 f$2 f$3 f$4))
define(backQuad, (b$1 b$4 b$3 b$2))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters

convertToMeters 1;

define(r, 0.000625)

define(maxRefLevel, 2)
define(quot, calc(2**maxRefLevel) )
define(n, calc(16/quot))




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//Derived parameters

define(dx, calc(r/n))
define(dxMax, calc(r/16))

define(h, calc(17*r)) 
define(l,  calc(r))
define(b,  dx)
define(ny, calc(h/dx))
define(nx, calc(r/dx))
define(nz,  1)



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Parametric description

  vertices
  (
      (0.0  0.0  0.0) vlabel(p0)
      (r    0.0  0.0) vlabel(p1)
      (r    h   0.0) vlabel(p2)
      (0.0  h   0.0) vlabel(p3)

      (0.0  0.0  b) vlabel(p4)
      (r    0.0  b) vlabel(p5)
      (r    h   b) vlabel(p6)
      (0.0  h   b) vlabel(p7)

  );

  blocks 
  (
      hex (p0 p1 p2 p3 p4 p5 p6 p7) (nx ny nz) simpleGrading (1 1 1)

  );
    
  edges ( );
    
  patches 
  (
	symmetry left 
  ( 	
		(p4 p7 p3 p0)
	)
 
	patch bottom
  (  
		(p1 p5 p4 p0)
	)

	patch right
  ( 	
		(p1 p2 p6 p5)
  )

	patch top
  ( 	
		 (p3 p7 p6 p2)
	)

  symmetry front
  (   
     (p0 p1 p2 p3)
  )

  symmetry back
  (   
     (p6 p5 p4 p7)
  )


  );
    
 mergePatchPairs ( );
