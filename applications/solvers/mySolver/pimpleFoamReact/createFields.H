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

//
// ############################################ create reaction dict ############################################
//

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

// Optional output of derived fields

const Switch writeCFields
(
    reactions.getOrDefault<Switch>
    (
        "writeCFields",
        false
    )
);

const Switch writeRYFields
(
    reactions.getOrDefault<Switch>
    (
        "writeRYFields",
        false
    )
);

const Switch writePEqFields
(
    reactions.getOrDefault<Switch>
    (
        "writePEqFields",
        false
    )
);

const Switch writeSFields
(
    reactions.getOrDefault<Switch>
    (
        "writeSFields",
        false
    )
);

const Switch writePEqTotField
(
    reactions.getOrDefault<Switch>
    (
        "writePEqTotField",
        false
    )
);

const Switch writeSTotField
(
    reactions.getOrDefault<Switch>
    (
        "writeSTotField",
        false
    )
);

//
// ############################################ create species mass-fraction field ############################################
//

wordList species(reactions.lookup("species"));
word inertSpecie(reactions.lookup("inertSpecie"));

PtrList<volScalarField> Y(species.size());
PtrList<volScalarField> RY(species.size());

forAll(species, i)
{
    Info<< "Reading species field Y_" << species[i] << endl;
    
    Y.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Y_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    RY.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RY_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writeRYFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );
}

//
// ############################################ convert to concentration ############################################
//

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
        "M_" + species[i],
        dimMass/dimMoles,
        molarMassDict.get<scalar>(species[i])
    );

    c.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "c_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writeCFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            rho*Yi/Mi
        )
    );
}

//
// ############################################ saturation diagnostics  ############################################
//

// Henry-law saturation diagnostics setup
dimensionedScalar pAbsRef
(
    "absolutePressureReference",
    dimPressure,
    reactions.getOrDefault<scalar>("absolutePressureReference", 101325.0)
);

volScalarField pEq_tot
(
    IOobject
    (
        "pEq_tot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writePEqTotField
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0.0)
);

volScalarField S_tot
(
    IOobject
    (
        "S_tot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeSTotField
            ? IOobject::AUTO_WRITE
            : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
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

dimensionedScalar pSmall
(
    "pSmall",
    dimPressure,
    SMALL
);

volScalarField pAbs
(
    IOobject
    (
        "pAbs_tmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pAbsRef + rho*p
);

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
                writePEqFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
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
                writeSFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            // pEqHenry[hi]/(pAbsRef + rho*p)
            pEqHenry[hi]/max(pAbs, pSmall)
        )
    );

    pEq_tot += pEqHenry[hi];
}
// S_tot = pEq_tot/(pAbsRef + rho*p);
S_tot = pEq_tot/max(pAbs, pSmall);

pEq_tot.correctBoundaryConditions();
S_tot.correctBoundaryConditions();


#include "createMRF.H"
#include "createFvOptions.H"

