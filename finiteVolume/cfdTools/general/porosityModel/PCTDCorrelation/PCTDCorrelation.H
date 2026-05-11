/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::porosityModels::PCTDCorrelation

Description
    Darcy-Forchheimer law porosity model, given by:

        \f[
            S = - (\mu d + \frac{\rho |U|}{2} f) U
        \f]

    where
    \vartable
        d        | Darcy coefficient [1/m2]
        f        | Forchheimer coefficient [1/m]
    \endvartable

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The orientation of the porous region is defined with the same notation as
    a co-ordinate system, but only a Cartesian co-ordinate system is valid.

SourceFiles
    PCTDCorrelation.C
    PCTDCorrelationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef PCTDCorrelation_H
#define PCTDCorrelation_H

#include "porosityModel.H"
#include "dimensionedTensor.H"

#include <vector>
#include <numbers>
#include "fvMatrix.H"
#include "vector.H"
#include "scalar.H"
#include <numeric> // For std::accumulate

#include <Eigen/Dense>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class PCTDCorrelation Declaration
\*---------------------------------------------------------------------------*/

class PCTDCorrelation
:
    public porosityModel
{
    // Private data

        //- Empirical L/D ratio
        scalar Nr_;
        //- Empirical Geometrical parameter
        scalar D_;
         //- Empirical F parameter
        scalar PD_;
         //- Empirical Re/U
        scalar Dw_;
        scalar H_;
        scalar W_;
        scalar r_;
		scalar L_active_;
		scalar L_porous_;

        scalar	m_ ;
        scalar    cbL1_  ;
        scalar    cbL2_  ;
        scalar    cbT1_  ;
        scalar    cbT2_  ;
        scalar    a_     ;
        scalar    b_     ;
        scalar    cwT1_  ;
        scalar    cwT2_  ;
        scalar    cwT3_  ;
        scalar    cwT4_  ;
        scalar    WmT_   ;
        scalar    cwL1_  ;
        scalar    cwL2_  ;
        scalar    WmL_   ;
        scalar    gamma_ ;
        scalar    lambda_;
	//- Name of density field
        word rhoName_;

        //- Name of dynamic viscosity field
        word muName_;

        //- Name of kinematic viscosity field
        word nuName_;

    // Private Member Functions

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        PCTDCorrelation(const PCTDCorrelation&);

        //- Disallow default bitwise assignment
        void operator=(const PCTDCorrelation&);


public:

    //- Runtime type information
    TypeName("PCTDCorrelation");

    //- Constructor
    PCTDCorrelation
    (
        const word& name,
    	const fvMesh& mesh,
    	const dictionary& dict,
    	const dictionary& coeffDict,
    	const word& cellZoneName
    );

    //- Destructor
    virtual ~PCTDCorrelation();


    // Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PCTDCorrelationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
