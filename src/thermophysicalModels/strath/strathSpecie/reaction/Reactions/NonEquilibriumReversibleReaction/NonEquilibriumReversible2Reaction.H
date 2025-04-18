/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::NonEquilibriumReversible2Reaction

Description
    Simple extension of Reaction2 to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    NonEquilibriumReversible2Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef NonEquilibriumReversible2Reaction_H
#define NonEquilibriumReversible2Reaction_H

#include "Reaction2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class NonEquilibriumReversible2Reaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class NonEquilibriumReversible2Reaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate fk_;
        ReactionRate rk_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const NonEquilibriumReversible2Reaction
                <ReactionType, ReactionThermo, ReactionRate>&
        );


public:

    //- Runtime type information
    TypeName("nonEquilibriumReversible");


    // Constructors

        //- Construct from components
        NonEquilibriumReversible2Reaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& forwardReactionRate,
            const ReactionRate& reverseReactionRate
        );

        //- Construct as copy given new speciesTable
        NonEquilibriumReversible2Reaction
        (
            const NonEquilibriumReversible2Reaction
                <ReactionType, ReactionThermo, ReactionRate>&,
            const speciesTable& species
        );

        //- Construct from Istream
        NonEquilibriumReversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        NonEquilibriumReversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<ReactionType<ReactionThermo> > clone() const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new NonEquilibriumReversible2Reaction
                    <ReactionType, ReactionThermo, ReactionRate>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<ReactionType<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new NonEquilibriumReversible2Reaction
                    <ReactionType, ReactionThermo, ReactionRate>
                    (*this, species)
            );
        }


    //- Destructor
    virtual ~NonEquilibriumReversible2Reaction()
    {}


    // Member Functions

        // NonEquilibriumReversible2Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NonEquilibriumReversible2Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
