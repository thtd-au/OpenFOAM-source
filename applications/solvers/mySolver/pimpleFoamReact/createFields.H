#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Reading reactions dictionary\n" << endl;

IOdictionary reactions
(
    IOobject
    (
        "reactions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

wordList species(reactions.lookup("species"));
word inertSpecie(reactions.lookup("inertSpecie"));

PtrList<volScalarField> Y(species.size());

forAll(species, i)
{
    Info<< "Reading species field " << species[i] << endl;

    Y.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                species[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

const dictionary& molarMassDict = reactions.subDict("molarMass");

PtrList<volScalarField> c(Y.size());

const dimensionedScalar rho
(
    "rho",
    dimDensity,
    laminarTransport
);

forAll(Y, i)
{
    const volScalarField& Yi = Y[i];

    dimensionedScalar Mi
    (
        "M_" + Yi.name(),
        dimMass/dimMoles,
        molarMassDict.get<scalar>(Yi.name())
    );

    c.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "c_" + Yi.name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            rho*Yi/Mi
        )
    );
}

// Henry-law saturation diagnostics setup
dimensionedScalar pAbsRef
(
    "absolutePressureReference",
    dimPressure,
    reactions.getOrDefault<scalar>("absolutePressureReference", 101325.0)
);

wordList henrySpecies(0);
labelList henrySpeciesIndex(0);
scalarList henryConstants(0);

if (reactions.found("henryConstants") || reactions.found("henryConstant"))
{
    const dictionary& henryDict = reactions.found("henryConstants")
      ? reactions.subDict("henryConstants")
      : reactions.subDict("henryConstant");

    henrySpecies = henryDict.toc();
    henrySpeciesIndex.setSize(henrySpecies.size(), -1);
    henryConstants.setSize(henrySpecies.size(), 0.0);

    Info<< "Creating Henry saturation diagnostics for species:" << nl;

    forAll(henrySpecies, hi)
    {
        const word& specieName = henrySpecies[hi];
        const label speciei = species.find(specieName);

        if (speciei < 0)
        {
            FatalIOErrorInFunction(reactions)
                << "Species " << specieName
                << " is listed in henryConstants but not in species"
                << exit(FatalIOError);
        }

        const scalar H = henryDict.get<scalar>(specieName);

        if (H <= SMALL)
        {
            FatalIOErrorInFunction(reactions)
                << "Henry constant for species " << specieName
                << " must be positive"
                << exit(FatalIOError);
        }

        henrySpeciesIndex[hi] = speciei;
        henryConstants[hi] = H;

        Info<< "  " << specieName << ": H = " << H
            << " mol/(m3 Pa)" << nl;
    }

    Info<< "absolutePressureReference = " << pAbsRef.value()
        << " Pa" << nl << endl;
}
else
{
    Info<< "Henry saturation diagnostics: none" << nl << endl;
}

PtrList<volScalarField> pEqHenry(henrySpecies.size());
PtrList<volScalarField> saturation(henrySpecies.size());

forAll(henrySpecies, hi)
{
    const label speciei = henrySpeciesIndex[hi];
    const word& specieName = henrySpecies[hi];

    dimensionedScalar Hi
    (
        "H_" + specieName,
        dimMoles/dimVolume/dimPressure,
        henryConstants[hi]
    );

    pEqHenry.set
    (
        hi,
        new volScalarField
        (
            IOobject
            (
                "pEq_" + specieName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            c[speciei]/Hi
        )
    );

    saturation.set
    (
        hi,
        new volScalarField
        (
            IOobject
            (
                "S_" + specieName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pEqHenry[hi]/(pAbsRef + rho*p)
        )
    );
}

#include "createMRF.H"
#include "createFvOptions.H"

