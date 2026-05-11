/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "PCTDCorrelation.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(PCTDCorrelation, 0);
        addToRunTimeSelectionTable(porosityModel, PCTDCorrelation, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::PCTDCorrelation::PCTDCorrelation
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& cellZoneName
)
:
    porosityModel(name, mesh, dict, coeffDict, cellZoneName),
    Nr_(coeffDict.lookup<scalar>("Nr")),
    D_(coeffDict.lookup<scalar>("D")),
    PD_(coeffDict.lookup<scalar>("PD")),    
    Dw_(coeffDict.lookup<scalar>("Dw")),
    H_(coeffDict.lookup<scalar>("H")),   
    W_(coeffDict.lookup<scalar>("W")),   
    r_(coeffDict.lookupOrDefault<scalar>("r", 0)), 
    L_active_(coeffDict.lookup<scalar>("L_active")), 
    L_porous_(coeffDict.lookup<scalar>("L_porous")), 


    m_     (coeffDict.lookupOrDefault<scalar>("m"      , 0.18)), 
    cbL1_  (coeffDict.lookupOrDefault<scalar>("cbL1"   , 700)),  
    cbL2_  (coeffDict.lookupOrDefault<scalar>("cbL2"   , 0)),    
    cbT1_  (coeffDict.lookupOrDefault<scalar>("cbT1"   , 10000)),
    cbT2_  (coeffDict.lookupOrDefault<scalar>("cbT2"   , 0)),    
    a_     (coeffDict.lookupOrDefault<scalar>("a"      , -6.9)), 
    b_     (coeffDict.lookupOrDefault<scalar>("b"      , 12.0)), 
    cwT1_  (coeffDict.lookupOrDefault<scalar>("cwT1"   , 15.2)), 
    cwT2_  (coeffDict.lookupOrDefault<scalar>("cwT2"   , -48.0)),
    cwT3_  (coeffDict.lookupOrDefault<scalar>("cwT3"   , 148.6)),
    cwT4_  (coeffDict.lookupOrDefault<scalar>("cwT4"   , -0.547)),
    WmT_   (coeffDict.lookupOrDefault<scalar>("WmT"    , 8.8)),  
    cwL1_  (coeffDict.lookupOrDefault<scalar>("cwL1"   , 0.80)), 
    cwL2_  (coeffDict.lookupOrDefault<scalar>("cwL2"   , 1.20)), 
    WmL_   (coeffDict.lookupOrDefault<scalar>("WmL"    , 0)),    
    gamma_ (coeffDict.lookupOrDefault<scalar>("gamma"  , 0.362)),
    lambda_(coeffDict.lookupOrDefault<scalar>("lambda" , 6.7)),        
    rhoName_(coeffDict.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffDict.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffDict.lookupOrDefault<word>("nu", "nu"))
{
    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::PCTDCorrelation::~PCTDCorrelation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::PCTDCorrelation::calcTransformModelData()
{
    

}


void Foam::porosityModels::PCTDCorrelation::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;

}


void Foam::porosityModels::PCTDCorrelation::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
        }
    }

}


void Foam::porosityModels::PCTDCorrelation::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);


}


void Foam::porosityModels::PCTDCorrelation::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }

}





// ************************************************************************* //
