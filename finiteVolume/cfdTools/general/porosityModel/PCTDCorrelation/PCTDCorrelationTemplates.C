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
void Foam::porosityModels::PCTDCorrelation::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{

        #include "geometryParameters"
        #include "modelParameters"


            const scalar RebL  = pow(10., (cbL2 * (PD-1)) )* cbL1;        // Laminar boundary
            const scalar RebT  = pow(10,(cbT2 * (PD-1))) * cbT1;        // Turbulent boundary  

        Info << "RebL : " << RebL << " | RebT : " << RebT << endl; 


        std::vector<scalar> aL_PDLow ={26.0, 26.18 , 26.98}; 
        std::vector<scalar> aL_PDHigh ={62.97 , 44.4, 87.26}; 

        std::vector<scalar> aT_PDLow  ={0.09378, 0.09377, 0.1004}; 
        std::vector<scalar> aT_PDHigh={0.1458, 0.1430, 0.1499}; 

        std::vector<scalar> bL_PDLow  ={888.2, 554.5, 1636.0}; 
        std::vector<scalar> bL_PDHigh={ 216.9, 256.7, 38.59}; 

        std::vector<scalar> bT_PDLow  ={1.398, 0.8732, 1.625}; 
        std::vector<scalar> bT_PDHigh={0.03632, 0.04199, 0.006706}; 

        std::vector<scalar> cL_PDLow  ={-3334.0 , -1480, -10050}; 
        std::vector<scalar> cL_PDHigh={-190.2, -267.6, -55.12}; 

        std::vector<scalar> cT_PDLow  ={-8.664, -3.341, -11.85}; 
        std::vector<scalar> cT_PDHigh={-0.03333, -0.04428, -0.009567}; 

        std::vector<scalar> CfL_bare(3); 
        std::vector<scalar> CfT_bare(3); 

        // 1.- Interarior sub-channels

        if (PD<=1.1) {
            CfL_bare[0] = aL_PDLow[0] + bL_PDLow[0] * (PD-1) + cL_PDLow[0] * pow((PD-1),2);
            CfT_bare[0] = aT_PDLow[0] + bT_PDLow[0] * (PD-1) + cT_PDLow[0] * pow((PD-1),2);
        } else {
            CfL_bare[0] = aL_PDHigh[0] + bL_PDHigh[0] * (PD-1) + cL_PDHigh[0] * pow((PD-1),2);
            CfT_bare[0] = aT_PDHigh[0] + bT_PDHigh[0] * (PD-1) + cT_PDHigh[0] * pow((PD-1),2);

        }


        // 2. and 3.- Edge and corder sub-channels
        if (WD<=1.1) {
            CfL_bare[1] = aL_PDLow[1] + bL_PDLow[1] * (WD-1) + cL_PDLow[1] * pow((WD-1),2);
            CfT_bare[1] = aT_PDLow[1] + bT_PDLow[1] * (WD-1) + cT_PDLow[1] * pow((WD-1),2);
            
            CfL_bare[2] = aL_PDLow[2] + bL_PDLow[2] * (WD-1) + cL_PDLow[2] * pow((WD-1),2);
            CfT_bare[2] = aT_PDLow[2] + bT_PDLow[2] * (WD-1) + cT_PDLow[2] * pow((WD-1),2);

        } else {
            CfL_bare[1] = aL_PDHigh[1] + bL_PDHigh[1] * (WD-1) + cL_PDHigh[1] * pow((WD-1),2);
            CfT_bare[1] = aT_PDHigh[1] + bT_PDHigh[1] * (WD-1) + cT_PDHigh[1] * pow((WD-1),2);
            
            CfL_bare[2] = aL_PDHigh[2] + bL_PDHigh[2] * (WD-1) + cL_PDHigh[2] * pow((WD-1),2);
            CfT_bare[2] = aT_PDHigh[2] + bT_PDHigh[2] * (WD-1) + cT_PDHigh[2] * pow((WD-1),2);
        }

        scalar WdT = (cwT1 + cwT2 * (Dw/D) + cwT3 * pow((Dw/D),2)) * pow((HD),cwT4);

        scalar WsT = a * log10(HD) + b; 

        scalar WdL = cwL1 * WdT;
        scalar WsL = cwL2 * WsT;


 //       Info << "CfL_bare | CfT_bare : " << CfL_bare[1] << " | " <<CfL_bare[2] << " || " << CfT_bare[1] << " | " <<CfT_bare[2] << endl; 
 //       Info << "WdT | WsT | WdL | WsL : " << WdT << " | " << WsT << " | " << WdL << " | " << WsL << endl; 
 //       Info << "W/D : " << WD << endl; 
        const labelList& cells = mesh_.cellZones()[zoneName_];

        forAll(cells, i)
        {
            const scalar Re =  rho[cells[i]]*(mag(U[cells[i]]))*De_bdl/gamma_bdl/mu[cells[i]]+VSMALL;
            scalar rhoCell = rho[cells[i]];
            scalar muCell = mu[cells[i]];

     //     Info << "Re: " <<Re <<endl; 
      
	          tensor C2=I;
	          #include "PCTD_turbulent"   
	          #include "PCTD_laminar"  
            if (Re>RebT)
            {
               
                const scalar CfbT = pow((1/sumSX_T),(2-m));
                //Info << "CfbT : " << CfbT << " | " << sumSX_T << endl;
                const scalar fbT=CfbT/pow(Re, m);
    //            Info << "fbT : " << fbT << endl;
                scalar LLDD = L_active/De_bdl/pow(gamma_bdl,2)/L_porous;
                     
                C2=fbT*LLDD*0.5*I ; //1 has to be changed in L
    //            Info << "LLDD : " << LLDD << endl;
    //            Info << "Re : " << Re <<endl;
    //            Info << "fbT : " << fbT <<endl;
	        }
	        else if (Re<RebL)
	        {



                 scalar CfbL = pow((1/sumSX_L),(2-m));
                 scalar fbL=CfbL/pow(Re, m);
     //            Info << "fbL : " << fbL <<endl;
                scalar LLDD = L_active/De_bdl/pow(gamma_bdl,2)/L_porous;
                          
                C2=fbL*LLDD*0.5*I ; //1 has to be changed in L
    //            Info << "LLDD : " << LLDD << endl;
    //            Info << "Re : " << Re <<endl;
    //            Info << "fbL : " << fbL <<endl;


	        }
	        else { 
	        	          #include "PCTD_transition"   
	        	          
	        	 scalar fbTransition = pow((1/sumSX_transition),(2));
    //            Info << "fbTransition : " << fbTransition <<endl;
                scalar LLDD = L_active/De_bdl/pow(gamma_bdl,2)/L_porous;
                          
                C2=fbTransition*LLDD*0.5*I ; //1 has to be changed in L
	        }
	    
	    const tensor dragCoeff = (rho[cells[i]]*mag(U[cells[i]]))*C2;
            
            const scalar isoDragCoeff = tr(dragCoeff);

            Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
	    Usource[cells[i]] -=
                V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]);

        
    }

}


template<class RhoFieldType>
void Foam::porosityModels::PCTDCorrelation::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{

        #include "geometryParameters"

        #include "modelParameters"

            const scalar RebL  = pow(10., (cbL2 * (PD-1)) )* cbL1;        // Laminar boundary
            const scalar RebT  = pow(10,(cbT2 * (PD-1))) * cbT1;        // Turbulent boundary  

        const labelList& cells = mesh_.cellZones()[zoneName_];

        forAll(cells, i)
        {
	    const scalar Re = De_bdl/gamma_bdl * (mag(U[cells[i]])+ VSMALL);
	    tensor C2=I;
	    
	    if (Re>1)
	    {
	      const scalar fact = 0; //GE*(64./(Re)*sqrt(FF)+0.0816/(pow(Re, 0.133))*pow(FF, 0.9335));
	      C2=0; //fact*LD*0.5*I;
	    }
	    else
	    {
	      C2=I;
	    }

            AU[cells[i]] += (rho[cells[i]]*mag(U[cells[i]]))*C2;
        }
 

}


// ************************************************************************* //
