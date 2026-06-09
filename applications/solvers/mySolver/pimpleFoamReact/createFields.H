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

#include "createMRF.H"
#include "createFvOptions.H"

