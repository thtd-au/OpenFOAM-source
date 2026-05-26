#ifndef faradaicMoleFractionFluxFvPatchScalarField_H
#define faradaicMoleFractionFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class faradaicMoleFractionFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    scalar stoichCoeff_;
    scalar nElectrons_;
    scalar currentDensity_;

public:

    TypeName("faradaicMoleFractionFlux");       // name used in type field for BC in 0/species 

    faradaicMoleFractionFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    faradaicMoleFractionFluxFvPatchScalarField
    (
        const faradaicMoleFractionFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    faradaicMoleFractionFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    faradaicMoleFractionFluxFvPatchScalarField
    (
        const faradaicMoleFractionFluxFvPatchScalarField&
    );

    faradaicMoleFractionFluxFvPatchScalarField
    (
        const faradaicMoleFractionFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMoleFractionFluxFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMoleFractionFluxFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

