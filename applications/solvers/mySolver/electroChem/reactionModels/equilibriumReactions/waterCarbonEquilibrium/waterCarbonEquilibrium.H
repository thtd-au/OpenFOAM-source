#ifndef waterCarbonEquilibrium_H
#define waterCarbonEquilibrium_H

#include "fvMesh.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam
{

class waterCarbonEquilibrium
{
private:

    word name_;
    fvMesh& mesh_;

    word HplusName_;
    word OHminusName_;
    word CO2Name_;
    word CarbonateName_;
    word BicarbonateName_;

    // Chemistry parameters. Concentrations use mol/m3 and time uses s.
    scalar Kw_;
    scalar K2_;
    scalar K12_;
    scalar rf1_;
    scalar rb1_;
    scalar rf3_;
    scalar rb3_;

    // Local Newton controls read from system/fvSolution.
    label maxIter_;
    scalar absTol_;
    scalar relTol_;

    volScalarField& Hplus_;
    volScalarField& OHminus_;
    volScalarField& CO2_;
    volScalarField& Carbonate_;
    volScalarField& Bicarbonate_;

    tmp<volScalarField> zeroSource
    (
        const word& speciesName
    ) const;

public:

    waterCarbonEquilibrium
    (
        const word& name,
        fvMesh& mesh,
        const dictionary& reactionDict,
        const dictionary& solverControls
    );

    // Slow carbon pathways r1 and r3 enter the transported kinetic source.
    tmp<volScalarField> source
    (
        const word& speciesName
    ) const;

    // Fast rW/r2 reaction-extent projection using a cell-local 2x2 Newton solve.
    void apply();

    const word& name() const
    {
        return name_;
    }

    const word& HplusName() const
    {
        return HplusName_;
    }

    const word& OHminusName() const
    {
        return OHminusName_;
    }
};

} // End namespace Foam

#endif

