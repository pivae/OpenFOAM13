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
    Foam::porosityModels::RhemeCorrelation

Description
    Rheme correlation law porosity model, given by:

        \f[
            S = - ( \frac{\rho |U|}{2} f L/Dh) U 
        \f]
        
    with f being:
    
        
        \f[
        
            f=GE [-64/Re sqrt(FF) + 0.0816 (FF^0.93305)/Re^0.133]
        \f]
    where: 
        \f[        
            GE=N_rods \pi (D+D_w)/S_t with S_t being the total wetted perimeter in the FA
        \f]
        
        \f[
        
            FF=(P/D)^0.5 + [7.6 d_m/H (P/D)^2]^2.16
        \f]


    The input parameters are:
    \vartable
        LD          | L/Dh/porosity^2 |The porosity is used to take into account the actual bulk velocity in the FA
        GE          | S_b/S_t
        FF          | Geometry factor
        RePerU      | Deq/porosity
    \endvartable

    Based on: 
    Rehme, K. (1973). Pressure drop correlations for fuel element spacers. 
    Nuclear technology, 17(1), 15-23.
    https://doi.org/10.13182/NT73-A31250

    The orientation of the porous region is defined with the same notation as
    a co-ordinate system, but only a Cartesian co-ordinate system is valid.

SourceFiles
    RhemeCorrelation.C
    RhemeCorrelationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef RhemeCorrelation_H
#define RhemeCorrelation_H

#include "porosityModel.H"
#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class RhemeCorrelation Declaration
\*---------------------------------------------------------------------------*/

class RhemeCorrelation
:
    public porosityModel
{
    // Private data

        //- Empirical L/D ratio
        scalar LD_;

        //- Empirical Geometrical parameter
        scalar GE_;

         //- Empirical F parameter
        scalar FF_;
	
         //- Empirical Re/U
        scalar RePerU_;

		
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
        RhemeCorrelation(const RhemeCorrelation&);

        //- Disallow default bitwise assignment
        void operator=(const RhemeCorrelation&);


public:

    //- Runtime type information
    TypeName("RhemeCorrelation");

    //- Constructor
    RhemeCorrelation
    (
        const word& name,
    	const fvMesh& mesh,
    	const dictionary& dict,
    	const dictionary& coeffDict,
    	const word& cellZoneName
    );

    //- Destructor
    virtual ~RhemeCorrelation();


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
    #include "RhemeCorrelationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
