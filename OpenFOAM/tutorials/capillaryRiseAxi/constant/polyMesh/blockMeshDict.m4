/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


changecom(//)changequote([,])
//define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')])
define(calc, [esyscmd(perl -e 'print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(int_, [esyscmd(perl -e 'print ($1)')])

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

define(half_angle, calc(2.5*3.14159265/180.0))

define(dx, calc(r/n))
define(dxMax, calc(r/16))

define(h, calc(17*r) )
define(l,  calc(r) )
define(b,  calc(r * sin(half_angle)))
define(rx, calc(r * cos(half_angle)))
define(ny, calc(h/dx))
define(nx, calc(r/dx))
define(nz,  1)



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Parametric description

  vertices
  (
      (0.0  0.0  0.0) vlabel(p0)
      (rx   0.0 -b) vlabel(p1)
      (rx   h   -b) vlabel(p2)
      (0.0  h    0.0) vlabel(p3)

      (0.0  0.0  0.0) vlabel(p4)
      (rx   0.0  b) vlabel(p5)
      (rx   h    b) vlabel(p6)
      (0.0  h    0.0) vlabel(p7)

  );

  blocks
  (
      hex (p0 p1 p2 p3 p4 p5 p6 p7) (nx ny nz) simpleGrading (1 1 1)

  );

  edges ( );

boundary
(
    left
    {
        type symmetry;
        faces
        (
                (p4 p7 p3 p0)
        );
    }

    bottom
    {
        type patch;
        faces
        (
                (p1 p5 p4 p0)
        );
    }

    right
    {
        type patch;
        faces
        (
            (p1 p2 p6 p5)
        );
    }

    top
    {
        type patch;
        faces
        (
            (p3 p7 p6 p2)
        );
    }

    front
    {
        type wedge;
        faces
        (
           (p0 p1 p2 p3)
        );
    }

    back
    {
        type wedge;
        faces
        (
           (p6 p5 p4 p7)
        );
    }

  );

 mergePatchPairs ( );
