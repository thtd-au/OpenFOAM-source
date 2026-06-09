Info<< "Updating Henry saturation diagnostics" << endl;

if (henrySpecies.size())
{
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

        pEqHenry[hi] = c[speciei]/Hi;
        pEqHenry[hi].correctBoundaryConditions();

        saturation[hi] = pEqHenry[hi]/max(pAbs, pSmall);
        saturation[hi].correctBoundaryConditions();

        // Info<< pEqHenry[hi].name()
        //     << " min/max = "
        //     << min(pEqHenry[hi]).value() << " / "
        //     << max(pEqHenry[hi]).value()
        //     << " Pa" << endl;

        Info<< saturation[hi].name()
            << " min/max = "
            << min(saturation[hi]).value() << " / "
            << max(saturation[hi]).value()
            << endl;
    }
}

