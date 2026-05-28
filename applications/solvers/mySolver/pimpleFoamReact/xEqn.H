volScalarField xt
(
    IOobject
    (
        "xt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// diagnostics 
volScalarField xTotal
(
    IOobject
    (
        "xTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

label inertIndex = species.find(inertSpecie);

forAll(x, i)
{
    if (i != inertIndex)
    {
        volScalarField& xi = x[i];

        dimensionedScalar Di
        (
            "D_" + xi.name(),
            dimArea/dimTime,
            diffusivityDict.get<scalar>(xi.name())
        );

        fvScalarMatrix xiEqn
        (
            fvm::ddt(xi)
          + fvm::div(phi, xi)
          - fvm::laplacian(Di, xi)
         ==
            fvOptions(xi)
        );

        xiEqn.relax();
        fvOptions.constrain(xiEqn);
        xiEqn.solve(xi.name());
        fvOptions.correct(xi);

        xi.max(0.0);
        // xi.min(1.0); // can hide conservation errors and is therefore not active 
        xi.correctBoundaryConditions();

        xt += xi;
    }
}

if (inertIndex >= 0)
{
    x[inertIndex] = scalar(1) - xt;
    x[inertIndex].max(0.0);
    x[inertIndex].correctBoundaryConditions();
}

forAll(x, i)
{
    Info<< x[i].name()
        << " min/max = "
        << min(x[i]).value() << " / "
        << max(x[i]).value() << endl;
    
        xTotal += x[i];
}

Info<< "sum non-inert mole fractions xt min/max = "
    << min(xt).value() << " / "
    << max(xt).value() << endl;

if (max(xt).value() > 1.0)
{
    WarningInFunction
        << "Sum of non-inert mole fractions exceeds 1. "
        << "Inert species will become negative before clipping."
        << endl;
}

Info<< "total mole fraction xTotal min/max = "
    << min(xTotal).value() << " / "
    << max(xTotal).value() << endl;

Info<< "mole fraction closure error |xTotal - 1| max = "
    << max(mag(xTotal - scalar(1))).value() << endl;
