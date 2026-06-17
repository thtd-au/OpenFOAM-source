Info<< "Updating species concentrations c_i = rho*Y_i/M_i" << endl;

const dictionary& molarMassDict = reactions.subDict("molarMass");

forAll(Y, i)
{
    const volScalarField& Yi = Y[i];

    dimensionedScalar Mi
    (
        "M_" + species[i],
        dimMass/dimMoles,
        molarMassDict.get<scalar>(species[i])
    );

    c[i] = rho*Yi/Mi;
    c[i].correctBoundaryConditions();

    Info<< c[i].name()
        << " min/max = "
        << min(c[i]).value() << " / "
        << max(c[i]).value()
        << " mol/m3" << endl;
}

