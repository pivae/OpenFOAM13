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
#include "RhemeCorrelation.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(RhemeCorrelation, 0);
        addToRunTimeSelectionTable(porosityModel, RhemeCorrelation, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::RhemeCorrelation::RhemeCorrelation
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& cellZoneName
)
:
    porosityModel(name, mesh, dict, coeffDict, cellZoneName),
    LD_(coeffDict.lookupOrDefault<scalar>("LD", 0)),
    GE_(coeffDict.lookupOrDefault<scalar>("GE", 0)),
    FF_(coeffDict.lookupOrDefault<scalar>("FF", 0)),    
    RePerU_(coeffDict.lookupOrDefault<scalar>("RePerU", 1)),    
    rhoName_(coeffDict.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffDict.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffDict.lookupOrDefault<word>("nu", "nu"))
{
    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::RhemeCorrelation::~RhemeCorrelation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::RhemeCorrelation::calcTransformModelData()
{
    

}


void Foam::porosityModels::RhemeCorrelation::calcForce
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


void Foam::porosityModels::RhemeCorrelation::correct
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


void Foam::porosityModels::RhemeCorrelation::correct
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


void Foam::porosityModels::RhemeCorrelation::correct
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
