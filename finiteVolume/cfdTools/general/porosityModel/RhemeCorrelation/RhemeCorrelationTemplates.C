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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::RhemeCorrelation::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{

    
    
        const scalar LD = LD_;
        const scalar GE = GE_;
        const scalar FF = FF_;
        const scalar RePerU = RePerU_;
	
	
        const labelList& cells = mesh_.cellZones()[zoneName_];

        forAll(cells, i)
        {
            const scalar Re =  rho[cells[i]]*(mag(U[cells[i]]))*RePerU/mu[cells[i]]+VSMALL;
	    tensor C2=I;
	    
	    if (Re>1)
	    {
	      const scalar fact = GE*(64./(Re)*sqrt(FF)+0.0816/(pow(Re, 0.133))*pow(FF, 0.9335));
	      C2=fact*LD*0.5*I;
	    }
	    else
	    {
	      C2=I;
	    }
	    
	    const tensor dragCoeff = (rho[cells[i]]*mag(U[cells[i]]))*C2;
            
            const scalar isoDragCoeff = tr(dragCoeff);

            Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
	    Usource[cells[i]] -=
                V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]);

        
    }

}


template<class RhoFieldType>
void Foam::porosityModels::RhemeCorrelation::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{

    
        const scalar LD = LD_;
        const scalar GE = GE_;
        const scalar FF = FF_;
        const scalar RePerU = RePerU_;     

        const labelList& cells = mesh_.cellZones()[zoneName_];

        forAll(cells, i)
        {
	    const scalar Re = RePerU * (mag(U[cells[i]])+ VSMALL);
	    tensor C2=I;
	    
	    if (Re>1)
	    {
	      const scalar fact = GE*(64./(Re)*sqrt(FF)+0.0816/(pow(Re, 0.133))*pow(FF, 0.9335));
	      C2=fact*LD*0.5*I;
	    }
	    else
	    {
	      C2=I;
	    }

            AU[cells[i]] += (rho[cells[i]]*mag(U[cells[i]]))*C2;
        }
 

}


// ************************************************************************* //
