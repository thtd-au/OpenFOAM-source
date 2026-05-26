Info<< "Solving species equations" << endl;

volScalarField Yt
(
    IOobject
    (
        "Yt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

label inertIndex = species.find(inertSpecie);

forAll(Y, i)
{
    if (i != inertIndex)
    {
        volScalarField& Yi = Y[i];

        dimensionedScalar Di
        (
            "D_" + Yi.name(),
            dimArea/dimTime,
            diffusivityDict.get<scalar>(Yi.name())
        );

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
          - fvm::laplacian(Di, Yi)
         ==
            fvOptions(Yi)
        );

        YiEqn.relax();
        fvOptions.constrain(YiEqn);
        YiEqn.solve(Yi.name());
        fvOptions.correct(Yi);

        Info<< Yi.name()
        << " before clipping min/max = "
        << min(Yi).value() << " / "
        << max(Yi).value() << endl;

        Yi.max(0.0);
        // Yi.min(1.0); // can hide conservation errors and is therefore not active 

        Info<< Yi.name()
        << " after clipping min/max = "
        << min(Yi).value() << " / "
        << max(Yi).value() << endl;

        Yt += Yi;
    }
}

if (inertIndex >= 0)
{
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    Y[inertIndex].correctBoundaryConditions();
}

