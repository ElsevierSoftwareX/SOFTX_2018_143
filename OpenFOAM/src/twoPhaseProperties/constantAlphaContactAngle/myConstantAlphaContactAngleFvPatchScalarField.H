/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::myConstantAlphaContactAngleFvPatchScalarField

Description
    A constant alphaContactAngle scalar boundary condition.

See also
    Foam::alphaContactAngleFvPatchScalarField
    Foam::temperatureDependentAlphaContactAngleFvPatchScalarField

SourceFiles
    myConstantAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myConstantAlphaContactAngleFvPatchScalarField_H
#define myConstantAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class myConstantAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class myConstantAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;


public:

    //- Runtime type information
    TypeName("myConstantAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        myConstantAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myConstantAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myConstantAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        myConstantAlphaContactAngleFvPatchScalarField
        (
            const myConstantAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myConstantAlphaContactAngleFvPatchScalarField
        (
            const myConstantAlphaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myConstantAlphaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myConstantAlphaContactAngleFvPatchScalarField
        (
            const myConstantAlphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myConstantAlphaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
