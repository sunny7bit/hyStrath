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
    Foam::VTRelaxationModels::ShatalovMWP

Description
    Millikan-White-Park V-T Relaxation model with Shatalov correlation
    for O2-O2 interactions.
    See:
	"Investigation of oxygen dissociation and vibrational relaxation
	at temperatures 4000-10800 K".
	L. B. Ibraguimova, A. L. Sergievskaya, V. Yu. Levashov,
	O. P. Shatalov, Yu. V. Tunik, and I. E. Zabelinskii.
	http://dx.doi.org/10.1063/1.4813070
    
SourceFiles
    ShatalovMWP.C

\*---------------------------------------------------------------------------*/

#ifndef ShatalovMWP_H
#define ShatalovMWP_H

#include "VTRelaxationModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace VTRelaxationModels
{

/*---------------------------------------------------------------------------*\
                           Class ShatalovMWP Declaration
\*---------------------------------------------------------------------------*/

class ShatalovMWP
:
    public VTRelaxationModel
{
    // Private data

        // Coefficient
        scalar species1_, species2_;
        
        // Shatalov's correction or MWP VT energy exchange?
        bool ShatON_;
        
        // Coefficients
        scalar A12_, B12_; 
        
        // Molecular weight, species i [g/mol]
        scalar W1_;   
        
        // Offset term in the exponential
        scalar offset_;   
        
        // Coefficients related to Park's correction
        scalar sigma1_, sigma2_;
        
        // Temperature limits for Shatalov's empirical correction
	scalar Tlow_, Thigh_; 
        
        // Coefficients related to Shatalov's empirical correction
        FixedList<scalar, 4> empir_;


public:

    //- Runtime type information
    TypeName("ShatalovMWP");


    // Constructors

        //- construct from components
        ShatalovMWP
        (
            const word& name1,
            const word& name2,
            const label& lname1,
            const label& lname2,
            const dictionary& dict2T,
            const dictionary& dictThermoPhy,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD
        );


    // Destructor

        ~ShatalovMWP()
        {}


    // Member Functions

        //- V-T relaxation time [s]
        virtual tmp<volScalarField> tauVT() const;
         
        //- V-T relaxation time for patch [s]
        virtual tmp<scalarField> tauVT
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tt,
            const PtrList<scalarField>& Tv,
            const PtrList<scalarField>& nD 
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VTRelaxationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
