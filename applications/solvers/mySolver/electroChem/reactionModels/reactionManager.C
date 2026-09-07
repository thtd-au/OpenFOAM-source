#include "reactionManager.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

namespace Foam
{

reactionManager::reactionManager
(
    fvMesh& mesh,
    const dictionary& reactionsDict
)
:
    mesh_(mesh),
    reactionsDict_(reactionsDict),
    reversibleKineticReactions_(),
    monodReactions_(),
    waterEquilibriumReactions_(),
    waterCarbonEquilibriumReactions_()
{
    readKineticReactions();
    readEquilibriumReactions();

    Info<< "Reaction manager loaded:" << nl
        << "    reversible kineti   c reactions: "
        << reversibleKineticReactions_.size() << nl
        << "    multiplicative Monod reactions: "
        << monodReactions_.size() << nl
        << "    water equilibrium reactions: "
        << waterEquilibriumReactions_.size() << nl
        << "    water-carbon equilibrium reactions: "
        << waterCarbonEquilibriumReactions_.size() << nl;
}


void reactionManager::readKineticReactions()
{
    if (!reactionsDict_.found("kineticReactions"))
    {
        Info<< "No kineticReactions dictionary found." << nl;

        monodReactions_.setSize(0);
        reversibleKineticReactions_.setSize(0);
        return;
    }

    const dictionary& kineticDict =
        reactionsDict_.subDict("kineticReactions");

    const wordList reactionNames =
        kineticDict.toc();

    label nMonod = 0;
    label nReversible = 0;

    // First pass: count each concrete model type.
    forAll(reactionNames, reactioni)
    {
        const word& reactionName =
            reactionNames[reactioni];

        if (!kineticDict.isDict(reactionName))
        {
            FatalIOErrorInFunction(kineticDict)
                << "Entry " << reactionName
                << " in kineticReactions is not a dictionary."
                << exit(FatalIOError);
        }

        const dictionary& reactionDict =
            kineticDict.subDict(reactionName);

        const word modelType
        (
            reactionDict.lookup("type")
        );

        if (modelType == "multiplicativeMonod")
        {
            ++nMonod;
        }
        else if (modelType == "reversibleKineticReaction")
        {
            ++nReversible;
        }
        else
        {
            FatalIOErrorInFunction(reactionDict)
                << "Unsupported kinetic reaction type "
                << modelType << nl
                << "Supported types are:" << nl
                << "    multiplicativeMonod" << nl
                << "    reversibleKineticReaction"
                << exit(FatalIOError);
        }
    }

    monodReactions_.setSize(nMonod);
    reversibleKineticReactions_.setSize(nReversible);

    label monodI = 0;
    label reversibleI = 0;

    // Second pass: construct the reactions.
    forAll(reactionNames, reactioni)
    {
        const word& reactionName =
            reactionNames[reactioni];

        const dictionary& reactionDict =
            kineticDict.subDict(reactionName);

        const word modelType
        (
            reactionDict.lookup("type")
        );

        if (modelType == "multiplicativeMonod")
        {
            monodReactions_.set
            (
                monodI++,
                new multiplicativeMonod
                (
                    reactionName,
                    mesh_,
                    reactionDict
                )
            );
        }
        else if (modelType == "reversibleKineticReaction")
        {
            reversibleKineticReactions_.set
            (
                reversibleI++,
                new reversibleKineticReaction
                (
                    reactionName,
                    mesh_,
                    reactionDict
                )
            );
        }
    }
}


void reactionManager::readEquilibriumReactions()
{
    if (!reactionsDict_.found("equilibriumReactions"))
    {
        Info<< "No equilibriumReactions dictionary found." << nl;

        waterEquilibriumReactions_.setSize(0);
        waterCarbonEquilibriumReactions_.setSize(0);
        return;
    }

    const dictionary& equilibriumDict =
        reactionsDict_.subDict("equilibriumReactions");

    const wordList reactionNames = equilibriumDict.toc();

    label nWater = 0;
    label nWaterCarbon = 0;

    // First pass: validate and count concrete equilibrium model types.
    forAll(reactionNames, reactioni)
    {
        const word& reactionName = reactionNames[reactioni];

        if (!equilibriumDict.isDict(reactionName))
        {
            FatalIOErrorInFunction(equilibriumDict)
                << "Entry " << reactionName
                << " in equilibriumReactions is not a dictionary."
                << exit(FatalIOError);
        }

        const dictionary& reactionDict =
            equilibriumDict.subDict(reactionName);

        const word modelType(reactionDict.lookup("type"));

        if (modelType == "waterEquilibrium")
        {
            ++nWater;
        }
        else if (modelType == "waterCarbonEquilibrium")
        {
            ++nWaterCarbon;
        }
        else
        {
            FatalIOErrorInFunction(reactionDict)
                << "Unsupported equilibrium reaction type "
                << modelType << nl
                << "Supported types are:" << nl
                << "    waterEquilibrium" << nl
                << "    waterCarbonEquilibrium"
                << exit(FatalIOError);
        }
    }

    waterEquilibriumReactions_.setSize(nWater);
    waterCarbonEquilibriumReactions_.setSize(nWaterCarbon);

    const dictionary* waterCarbonControlsPtr = nullptr;

    if (nWaterCarbon)
    {
        const dictionary& solutionDict = mesh_.solutionDict();

        if (!solutionDict.found("equilibriumProjection"))
        {
            FatalIOErrorInFunction(solutionDict)
                << "waterCarbonEquilibrium requires an equilibriumProjection "
                << "dictionary in system/fvSolution."
                << exit(FatalIOError);
        }

        const dictionary& projectionDict =
            solutionDict.subDict("equilibriumProjection");

        if (!projectionDict.found("waterCarbonEquilibrium"))
        {
            FatalIOErrorInFunction(projectionDict)
                << "Missing waterCarbonEquilibrium controls. Required entries:"
                << nl
                << "    maxIter" << nl
                << "    absTol" << nl
                << "    relTol"
                << exit(FatalIOError);
        }

        waterCarbonControlsPtr =
            &projectionDict.subDict("waterCarbonEquilibrium");
    }

    label waterI = 0;
    label waterCarbonI = 0;

    forAll(reactionNames, reactioni)
    {
        const word& reactionName = reactionNames[reactioni];
        const dictionary& reactionDict =
            equilibriumDict.subDict(reactionName);
        const word modelType(reactionDict.lookup("type"));

        if (modelType == "waterEquilibrium")
        {
            waterEquilibriumReactions_.set
            (
                waterI++,
                new waterEquilibrium
                (
                    reactionName,
                    mesh_,
                    reactionDict
                )
            );
        }
        else if (modelType == "waterCarbonEquilibrium")
        {
            waterCarbonEquilibriumReactions_.set
            (
                waterCarbonI++,
                new waterCarbonEquilibrium
                (
                    reactionName,
                    mesh_,
                    reactionDict,
                    *waterCarbonControlsPtr
                )
            );
        }
    }
}


tmp<volScalarField> reactionManager::zeroSource
(
    const word& speciesName
) const
{
    const word concentrationFieldName("c_" + speciesName);

    if (!mesh_.foundObject<volScalarField>(concentrationFieldName))
    {
        FatalErrorInFunction
            << "Cannot construct reaction source for species "
            << speciesName << nl
            << "Expected concentration field "
            << concentrationFieldName
            << exit(FatalError);
    }

    const volScalarField& concentration =
        mesh_.lookupObject<volScalarField>(concentrationFieldName);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "kineticSource_" + speciesName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                concentration.dimensions()/dimTime,
                0.0
            )
        )
    );
}

tmp<volScalarField> reactionManager::kineticSource
(
    const word& speciesName
) const
{
    tmp<volScalarField> tCombinedSource = zeroSource(speciesName);
    volScalarField& combinedSource = tCombinedSource.ref();

    forAll(monodReactions_, reactioni)
    {
        combinedSource +=
            monodReactions_[reactioni].source(speciesName);
    }

    forAll(reversibleKineticReactions_, reactioni)
    {
        combinedSource +=
            reversibleKineticReactions_[reactioni].source(speciesName);
    }

    // waterCarbonEquilibrium owns slow r1/r3 chemistry internally but exposes
    // it here so it enters the same transported kinetic source Rc.
    forAll(waterCarbonEquilibriumReactions_, reactioni)
    {
        combinedSource +=
            waterCarbonEquilibriumReactions_[reactioni].source(speciesName);
    }

    return tCombinedSource;
}


void reactionManager::applyEquilibriumReactions()
{
    forAll(waterEquilibriumReactions_, reactioni)
    {
        waterEquilibriumReactions_[reactioni].apply();
    }

    forAll(waterCarbonEquilibriumReactions_, reactioni)
    {
        waterCarbonEquilibriumReactions_[reactioni].apply();
    }
}

} // End namespace Foam

