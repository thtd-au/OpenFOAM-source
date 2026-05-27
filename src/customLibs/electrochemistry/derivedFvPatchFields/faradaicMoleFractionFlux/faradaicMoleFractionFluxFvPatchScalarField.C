#include "faradaicMoleFractionFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(faradaicMoleFractionFluxFvPatchScalarField, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMoleFractionFluxFvPatchScalarField,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMoleFractionFluxFvPatchScalarField,
    patchMapper
);


faradaicMoleFractionFluxFvPatchScalarField::
faradaicMoleFractionFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    stoichCoeff_(0),
    nElectrons_(1),
    currentDensity_(0)
{}


faradaicMoleFractionFluxFvPatchScalarField::
faradaicMoleFractionFluxFvPatchScalarField
(
    const faradaicMoleFractionFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


faradaicMoleFractionFluxFvPatchScalarField::
faradaicMoleFractionFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    stoichCoeff_(dict.get<scalar>("stoichCoeff")),
    nElectrons_(dict.get<scalar>("nElectrons")),
    currentDensity_(dict.get<scalar>("currentDensity"))
{
    // The base fixedGradient BC normally reads a user-supplied gradient.
    // Here the gradient is calculated later in updateCoeffs().  
    // So dummy gradient is created initialy and will be updated to actual gradient 
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );

    gradient() = 0.0;

    if (nElectrons_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "nElectrons must be positive" << exit(FatalIOError);
    }
}


faradaicMoleFractionFluxFvPatchScalarField::
faradaicMoleFractionFluxFvPatchScalarField
(
    const faradaicMoleFractionFluxFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


faradaicMoleFractionFluxFvPatchScalarField::
faradaicMoleFractionFluxFvPatchScalarField
(
    const faradaicMoleFractionFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


void faradaicMoleFractionFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const word specieName = internalField().name();

    IOdictionary reactions
    (
        IOobject
        (
            "reactions",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar F =
        reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329);

    const scalar cSol =
        reactions.get<scalar>("totalMolarConcentration");

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");

    if (!diffusivityDict.found(specieName))
    {
        FatalIOErrorInFunction(reactions)
            << "No diffusivity entry for species " << specieName
            << " in constant/reactions/diffusivity" << exit(FatalIOError);
    }

    const scalar D = diffusivityDict.get<scalar>(specieName);

    if (D <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Diffusivity for species " << specieName
            << " must be positive" << exit(FatalIOError);
    }

    if (cSol <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "totalMolarConcentration must be positive" << exit(FatalIOError);
    }

    gradient() =
        stoichCoeff_*currentDensity_
       /(nElectrons_*F*D*cSol);
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void faradaicMoleFractionFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("stoichCoeff", stoichCoeff_);
    os.writeEntry("nElectrons", nElectrons_);
    os.writeEntry("currentDensity", currentDensity_);
    writeEntry("value", os);
}

}

