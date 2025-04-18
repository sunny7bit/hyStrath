/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::noHEEnergyTransfer

Description
    noHEEnergyTransfer disables heavy-particles-electron energy transfer.
  
SourceFiles
    noHEEnergyTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef noHEEnergyTransfer_H
#define noHEEnergyTransfer_H

#include "relaxationTimeModelHE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class noHEEnergyTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class noHEEnergyTransfer
:
    public relaxationTimeModelHE
{  

protected:
  
    // Protected data
    
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;


    // Protected member functions
    

public:

    //- Runtime type information
    TypeName("noHEEnergyTransfer");

    
    // Constructors

        //- Construct from components
        noHEEnergyTransfer
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~noHEEnergyTransfer()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }
        
        // Correct the V-V relaxation times and V-V energy transfers
        virtual void correct();
                
        //- Read transportProperties dictionary
        virtual bool read();
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "noHEEnergyTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
