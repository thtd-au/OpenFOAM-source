#include "waterCarbonEquilibrium.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

namespace Foam
{

waterCarbonEquilibrium::waterCarbonEquilibrium
(
    const word& name,
    fvMesh& mesh,
    const dictionary& reactionDict,
    const dictionary& solverControls
)
:
    name_(name),
    mesh_(mesh),
    HplusName_(reactionDict.get<word>("Hplus")),
    OHminusName_(reactionDict.get<word>("OHminus")),
    CO2Name_(reactionDict.get<word>("CO2")),
    CarbonateName_(reactionDict.get<word>("Carbonate")),
    BicarbonateName_(reactionDict.get<word>("Bicarbonate")),
    Kw_(reactionDict.get<scalar>("Kw")),
    K2_(reactionDict.get<scalar>("K2")),
    K12_(reactionDict.get<scalar>("K12")),
    rf1_(reactionDict.get<scalar>("rf1")),
    rb1_(reactionDict.get<scalar>("rb1")),
    rf3_(reactionDict.get<scalar>("rf3")),
    rb3_(reactionDict.get<scalar>("rb3")),
    maxIter_(solverControls.get<label>("maxIter")),
    absTol_(solverControls.get<scalar>("absTol")),
    relTol_(solverControls.get<scalar>("relTol")),
    Hplus_
    (
        mesh_.lookupObjectRef<volScalarField>("c_" + HplusName_)
    ),
    OHminus_
    (
        mesh_.lookupObjectRef<volScalarField>("c_" + OHminusName_)
    ),
    CO2_
    (
        mesh_.lookupObjectRef<volScalarField>("c_" + CO2Name_)
    ),
    Carbonate_
    (
        mesh_.lookupObjectRef<volScalarField>("c_" + CarbonateName_)
    ),
    Bicarbonate_
    (
        mesh_.lookupObjectRef<volScalarField>("c_" + BicarbonateName_)
    )
{
    wordList modelSpecies(5);
    modelSpecies[0] = HplusName_;
    modelSpecies[1] = OHminusName_;
    modelSpecies[2] = CO2Name_;
    modelSpecies[3] = CarbonateName_;
    modelSpecies[4] = BicarbonateName_;

    forAll(modelSpecies, i)
    {
        for (label j = i + 1; j < modelSpecies.size(); ++j)
        {
            if (modelSpecies[i] == modelSpecies[j])
            {
                FatalIOErrorInFunction(reactionDict)
                    << "Reaction " << name_
                    << ": all configured species must be distinct. "
                    << "Duplicate species: " << modelSpecies[i]
                    << exit(FatalIOError);
            }
        }
    }

    if
    (
        Hplus_.dimensions() != OHminus_.dimensions()
     || Hplus_.dimensions() != CO2_.dimensions()
     || Hplus_.dimensions() != Carbonate_.dimensions()
     || Hplus_.dimensions() != Bicarbonate_.dimensions()
    )
    {
        FatalIOErrorInFunction(reactionDict)
            << "Reaction " << name_
            << ": all concentration fields must have identical dimensions."
            << exit(FatalIOError);
    }

    if (Kw_ <= 0.0 || K2_ <= 0.0 || K12_ <= 0.0)
    {
        FatalIOErrorInFunction(reactionDict)
            << "Reaction " << name_
            << ": Kw, K2 and K12 must be greater than zero."
            << exit(FatalIOError);
    }

    if (rf1_ < 0.0 || rb1_ < 0.0 || rf3_ < 0.0 || rb3_ < 0.0)
    {
        FatalIOErrorInFunction(reactionDict)
            << "Reaction " << name_
            << ": rf1, rb1, rf3 and rb3 must be non-negative."
            << exit(FatalIOError);
    }

    if (maxIter_ <= 0)
    {
        FatalIOErrorInFunction(solverControls)
            << "waterCarbonEquilibrium/maxIter must be greater than zero."
            << exit(FatalIOError);
    }

    if (absTol_ <= 0.0)
    {
        FatalIOErrorInFunction(solverControls)
            << "waterCarbonEquilibrium/absTol must be greater than zero."
            << exit(FatalIOError);
    }

    if (relTol_ < 0.0)
    {
        FatalIOErrorInFunction(solverControls)
            << "waterCarbonEquilibrium/relTol must be non-negative."
            << exit(FatalIOError);
    }
}


tmp<volScalarField> waterCarbonEquilibrium::zeroSource
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
            << "Expected concentration field " << concentrationFieldName
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
                name_ + "_source_" + speciesName,
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


tmp<volScalarField> waterCarbonEquilibrium::source
(
    const word& speciesName
) const
{
    // Slow pathways, positive left-to-right:
    //
    // r1: CO2 + H2O <-> HCO3- + H+
    //     r1 = rf1*cCO2 - (rb1/K12)*cHCO3*cH
    //
    // r3: CO2 + OH- <-> HCO3-
    //     r3 = rf3*cCO2*cOH - rb3*cHCO3

    tmp<volScalarField> tSource = zeroSource(speciesName);
    volScalarField& reactionSource = tSource.ref();

    const dimensionSet concentrationDimensions = Hplus_.dimensions();

    const dimensionedScalar rf1
    (
        "rf1_" + name_,
        dimless/dimTime,
        rf1_
    );

    const dimensionedScalar rb1
    (
        "rb1_" + name_,
        dimless/dimTime,
        rb1_
    );

    const dimensionedScalar K12
    (
        "K12_" + name_,
        concentrationDimensions,
        K12_
    );

    dimensionSet secondOrderRateConstantDimensions(dimless/dimTime);
    secondOrderRateConstantDimensions /= concentrationDimensions;

    const dimensionedScalar rf3
    (
        "rf3_" + name_,
        secondOrderRateConstantDimensions,
        rf3_
    );

    const dimensionedScalar rb3
    (
        "rb3_" + name_,
        dimless/dimTime,
        rb3_
    );

    tmp<volScalarField> tr1 =
        rf1*CO2_ - (rb1/K12)*Bicarbonate_*Hplus_;

    tmp<volScalarField> tr3 =
        rf3*CO2_*OHminus_ - rb3*Bicarbonate_;

    const volScalarField& r1 = tr1();
    const volScalarField& r3 = tr3();

    if (speciesName == CO2Name_)
    {
        reactionSource -= r1;
        reactionSource -= r3;
    }
    else if (speciesName == HplusName_)
    {
        reactionSource += r1;
    }
    else if (speciesName == OHminusName_)
    {
        reactionSource -= r3;
    }
    else if (speciesName == BicarbonateName_)
    {
        reactionSource += r1;
        reactionSource += r3;
    }
    // Carbonate has no slow-path source in the reduced reaction system.

    return tSource;
}


void waterCarbonEquilibrium::apply()
{
    scalarField& cH = Hplus_.primitiveFieldRef();
    scalarField& cOH = OHminus_.primitiveFieldRef();
    scalarField& cCarbonate = Carbonate_.primitiveFieldRef();
    scalarField& cBicarbonate = Bicarbonate_.primitiveFieldRef();

    // Internal relative determinant tolerance. This detects cancellation in
    // D = J11*J22 - J12*J21 without imposing an absolute concentration scale.
    const scalar determinantRelTol = 100.0*SMALL;

    scalar localInitialScaledMax = 0.0;
    scalar localFinalScaledMax = 0.0;
    scalar localInitialRawMax = 0.0;
    scalar localFinalRawMax = 0.0;
    label localMaxIterations = 0;
    label localUnconvergedCells = 0;

    scalar localWorstFailedScaled = -GREAT;
    label localWorstFailedCell = -1;

    scalar worstH = 0.0;
    scalar worstOH = 0.0;
    scalar worstBic = 0.0;
    scalar worstCar = 0.0;

    scalar worstFW = 0.0;
    scalar worstF2 = 0.0;

    scalar worstSW = 0.0;
    scalar worstS2 = 0.0;

    scalar worstScaledFW = 0.0;
    scalar worstScaledF2 = 0.0;

    scalar worstDXiW = 0.0;
    scalar worstDXi2 = 0.0;

    scalar worstD = 0.0;
    scalar worstScaledDet = 0.0;

    scalar worstAlpha = 1.0;

    forAll(cH, celli)
    {
        // c* is the transported state. Reaction extents are corrections from
        // this state, so the natural initial guess is xiW = xi2 = 0.
        const scalar h0 = cH[celli];
        const scalar oh0 = cOH[celli];
        const scalar bic0 = cBicarbonate[celli];
        const scalar car0 = cCarbonate[celli];

        if
        (
            h0 < 0.0
        || oh0 < 0.0
        || bic0 < 0.0
        || car0 < 0.0
        )
        {
            FatalErrorInFunction
                << "Reaction " << name_
                << " received a negative transported concentration "
                << "before equilibrium projection." << nl
                << "cell = " << celli << nl
                << "Hplus = " << h0 << " mol/m3" << nl
                << "OHminus = " << oh0 << " mol/m3" << nl
                << "Bicarbonate = " << bic0 << " mol/m3" << nl
                << "Carbonate = " << car0 << " mol/m3"
                << exit(FatalError);
        }

        scalar xiW = 0.0;
        scalar xi2 = 0.0;

        scalar initialRawResidual = 0.0;
        scalar initialScaledResidual = 0.0;
        scalar finalRawResidual = GREAT;
        scalar finalScaledResidual = GREAT;

        bool converged = false;
        label nIter = 0;

        scalar lastDXiW = 0.0;
        scalar lastDXi2 = 0.0;
        scalar lastD = 0.0;
        scalar lastScaledDet = 0.0;
        scalar lastAlpha = 1.0;

        while (true)
        {
            // Current Newton iterate concentrations c = c* + N*xi.
            const scalar h = h0 + xiW + xi2;
            const scalar oh = oh0 + xiW;
            const scalar bic = bic0 - xi2;
            const scalar car = car0 + xi2;

            // Nonlinear equilibrium constraints.
            const scalar FW = h*oh - Kw_;
            const scalar F2 = car*h - K2_*bic;

            const scalar rawResidual = max(mag(FW), mag(F2));

            const scalar SW = max(mag(h*oh), mag(Kw_));

            // absTol_ is also the internal scaling floor for the carbon
            // equation. It has the same dimensions as F2 when interpreted
            // in the case concentration units.
            const scalar S2 =
                max
                (
                    max(mag(car*h), mag(K2_*bic)),
                    absTol_
                );

            // const scalar scaledResidual = max(mag(FW)/SW, mag(F2)/S2);
            const scalar scaledFW = mag(FW)/SW;
            const scalar scaledF2 = mag(F2)/S2;
            const scalar scaledResidual = max(scaledFW, scaledF2);

            if (nIter == 0)
            {
                initialRawResidual = rawResidual;
                initialScaledResidual = scaledResidual;
            }

            finalRawResidual = rawResidual;
            finalScaledResidual = scaledResidual;

            if
            (
                rawResidual <= absTol_
             || scaledResidual <= relTol_
            )
            {
                converged = true;
                break;
            }

            if (nIter >= maxIter_)
            {
                break;
            }

            // Analytic Jacobian of [FW, F2] with respect to [xiW, xi2].
            const scalar J11 = h + oh;
            const scalar J12 = oh;
            const scalar J21 = car;
            const scalar J22 = car + h + K2_;

            const scalar determinantTerm1 = J11*J22;
            const scalar determinantTerm2 = J12*J21;
            const scalar D = determinantTerm1 - determinantTerm2;

            const scalar determinantScale =
                max(mag(determinantTerm1), mag(determinantTerm2));

            // If both determinant terms are exactly zero, D is singular.
            // Otherwise use a relative cancellation measure.
            if
            (
                determinantScale <= VSMALL
             || mag(D)/(determinantScale + VSMALL) <= determinantRelTol
            )
            {
                FatalErrorInFunction
                    << "Reaction " << name_
                    << " encountered a singular/near-singular Newton Jacobian."
                    << nl
                    << "cell = " << celli << nl
                    << "iteration = " << nIter << nl
                    << "xiW = " << xiW << " mol/m3" << nl
                    << "xi2 = " << xi2 << " mol/m3" << nl
                    << "Hplus = " << h << " mol/m3" << nl
                    << "OHminus = " << oh << " mol/m3" << nl
                    << "Bicarbonate = " << bic << " mol/m3" << nl
                    << "Carbonate = " << car << " mol/m3" << nl
                    << "FW = " << FW << " (mol/m3)^2" << nl
                    << "F2 = " << F2 << " (mol/m3)^2" << nl
                    << "J11 = " << J11 << " mol/m3" << nl
                    << "J12 = " << J12 << " mol/m3" << nl
                    << "J21 = " << J21 << " mol/m3" << nl
                    << "J22 = " << J22 << " mol/m3" << nl
                    << "D = " << D << " (mol/m3)^2" << nl
                    << "scaled determinant = "
                    << mag(D)/(determinantScale + VSMALL)
                    << exit(FatalError);
            }

            // Direct solution of J*dXi = -F for the 2x2 system.
            // const scalar dXiW = (-FW*J22 + J12*F2)/D;
            // const scalar dXi2 = (J21*FW - J11*F2)/D;

            const scalar dXiW = (-FW*J22 + J12*F2)/D;
            const scalar dXi2 = (J21*FW - J11*F2)/D;

            // ------------------------------------------------------------------ //
            // Positivity-preserving fraction-to-boundary step limiter
            //
            // Concentration changes associated with the Newton extent update:
            //
            //   dH   = dXiW + dXi2
            //   dOH  = dXiW
            //   dBic = -dXi2
            //   dCar = dXi2
            //
            // Reduce alpha only when the full Newton step would make one of the
            // equilibrium species negative.
            // ------------------------------------------------------------------ //

            const scalar dH   = dXiW + dXi2;
            const scalar dOH  = dXiW;
            const scalar dBic = -dXi2;
            const scalar dCar = dXi2;

            scalar alpha = 1.0;

            // Stay slightly inside the admissible region instead of landing
            // directly on c = 0.
            const scalar fractionToBoundary = 0.95;

            if (dH < 0.0)
            {
                // alpha = min(alpha,fractionToBoundary*h/(-dH));
                alpha = min(alpha,fractionToBoundary*max(h, scalar(0.0))/(-dH));
            }

            if (dOH < 0.0)
            {
                alpha = min(alpha,fractionToBoundary*oh/(-dOH));
            }

            if (dBic < 0.0)
            {
                alpha = min(alpha,fractionToBoundary*bic/(-dBic));
            }

            if (dCar < 0.0)
            {
                alpha = min(alpha,fractionToBoundary*car/(-dCar));
            }

            // alpha should never exceed the full Newton step.
            alpha = min(alpha, scalar(1.0));
            
            if (alpha <= VSMALL)
            {
                FatalErrorInFunction
                    << "Reaction " << name_
                    << " Newton step cannot remain inside the "
                    << "non-negative concentration region." << nl
                    << "cell = " << celli << nl
                    << "iteration = " << nIter << nl
                    << "alpha = " << alpha << nl
                    << "Hplus = " << h << nl
                    << "OHminus = " << oh << nl
                    << "Bicarbonate = " << bic << nl
                    << "Carbonate = " << car << nl
                    << "dXiW = " << dXiW << nl
                    << "dXi2 = " << dXi2
                    << exit(FatalError);
            }

            lastDXiW = dXiW;
            lastDXi2 = dXi2;
            lastD = D;
            lastScaledDet = mag(D)/(determinantScale + VSMALL);
            lastAlpha = alpha;

            xiW += alpha*dXiW;
            xi2 += alpha*dXi2;
            ++nIter;
        }

        localInitialRawMax = max(localInitialRawMax, initialRawResidual);
        localInitialScaledMax = max(localInitialScaledMax, initialScaledResidual);
        localFinalRawMax = max(localFinalRawMax, finalRawResidual);
        localFinalScaledMax = max(localFinalScaledMax, finalScaledResidual);
        localMaxIterations = max(localMaxIterations, nIter);

        // if (!converged)
        // {
        //     ++localUnconvergedCells;
        //     continue;
        // }
        if (!converged)
        {
            ++localUnconvergedCells;

            // Reconstruct the final unconverged Newton state
            const scalar h = h0 + xiW + xi2;
            const scalar oh = oh0 + xiW;
            const scalar bic = bic0 - xi2;
            const scalar car = car0 + xi2;

            const scalar FW = h*oh - Kw_;
            const scalar F2 = car*h - K2_*bic;

            const scalar SW = max(mag(h*oh), mag(Kw_));
            const scalar S2 = max(max(mag(car*h), mag(K2_*bic)),absTol_);

            const scalar scaledFW = mag(FW)/SW;
            const scalar scaledF2 = mag(F2)/S2;

            const scalar scaledResidual = max(scaledFW, scaledF2);

            if (scaledResidual > localWorstFailedScaled)
            {
                localWorstFailedScaled = scaledResidual;
                localWorstFailedCell = celli;

                worstH = h;
                worstOH = oh;
                worstBic = bic;
                worstCar = car;

                worstFW = FW;
                worstF2 = F2;

                worstSW = SW;
                worstS2 = S2;

                worstScaledFW = scaledFW;
                worstScaledF2 = scaledF2;

                worstDXiW = lastDXiW;
                worstDXi2 = lastDXi2;

                worstD = lastD;
                worstScaledDet = lastScaledDet;

                worstAlpha = lastAlpha;
            }

            continue;
        }

        // Only commit the projected state after the local Newton solve has
        // satisfied the convergence criterion.
        cH[celli] = h0 + xiW + xi2;
        cOH[celli] = oh0 + xiW;
        cBicarbonate[celli] = bic0 - xi2;
        cCarbonate[celli] = car0 + xi2;
    }

    scalar initialScaledMax = localInitialScaledMax;
    scalar finalScaledMax = localFinalScaledMax;
    scalar initialRawMax = localInitialRawMax;
    scalar finalRawMax = localFinalRawMax;
    label maxIterations = localMaxIterations;
    label unconvergedCells = localUnconvergedCells;

    reduce(initialScaledMax, maxOp<scalar>());
    reduce(finalScaledMax, maxOp<scalar>());
    reduce(initialRawMax, maxOp<scalar>());
    reduce(finalRawMax, maxOp<scalar>());
    reduce(maxIterations, maxOp<label>());
    reduce(unconvergedCells, sumOp<label>());

    Info<< name_
        << ": Solving for equilibrium"
        << ", initial scaled max = " << initialScaledMax
        << ", final scaled max = " << finalScaledMax
        << ", initial raw max = " << initialRawMax
        << ", final raw max = " << finalRawMax
        << ", max iterations = " << maxIterations
        << nl;

    if (unconvergedCells)
    {
        Info<< ", Unconverged cells = " << unconvergedCells;
    }

    Info<< nl;

    if (unconvergedCells > 0 && localWorstFailedCell >= 0)
    {
        Info<< nl
            << "Worst unconverged waterCarbonEquilibrium cell:" << nl
            << "    cell              = " << localWorstFailedCell << nl
            << "    Hplus             = " << worstH << nl
            << "    OHminus           = " << worstOH << nl
            << "    Bicarbonate       = " << worstBic << nl
            << "    Carbonate         = " << worstCar << nl
            << "    FW                 = " << worstFW << nl
            << "    F2                 = " << worstF2 << nl
            << "    SW                 = " << worstSW << nl
            << "    S2                 = " << worstS2 << nl
            << "    scaled FW          = " << worstScaledFW << nl
            << "    scaled F2          = " << worstScaledF2 << nl
            << "    last dXiW          = " << worstDXiW << nl
            << "    last dXi2          = " << worstDXi2 << nl
            << "    determinant        = " << worstD << nl
            << "    scaled determinant = " << worstScaledDet << nl
            << "    last alpha         = " << worstAlpha << nl
            << endl;
    }

    if (unconvergedCells)
    {
        FatalErrorInFunction
            << "Reaction " << name_
            << " Newton projection failed to converge in "
            << unconvergedCells << " cell(s)." << nl
            << "maxIter = " << maxIter_ << nl
            << "absTol = " << absTol_ << " (mol/m3)^2" << nl
            << "relTol = " << relTol_ << nl
            << "worst final raw residual = " << finalRawMax
            << " (mol/m3)^2" << nl
            << "worst final scaled residual = " << finalScaledMax
            << exit(FatalError);
    }

    // Internal cells are projected; configured patch conditions remain
    // authoritative, matching the existing equilibrium-model behaviour.
    Hplus_.correctBoundaryConditions();
    OHminus_.correctBoundaryConditions();
    Bicarbonate_.correctBoundaryConditions();
    Carbonate_.correctBoundaryConditions();

    // CO2 is intentionally not modified by the fast projection. It evolves
    // through transport and the slow r1/r3 kinetic source terms only.
}

} // End namespace Foam

