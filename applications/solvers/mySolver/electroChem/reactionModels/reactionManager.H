#ifndef reactionManager_H
#define reactionManager_H

#include "fvMesh.H"
#include "dictionary.H"
#include "PtrList.H"
#include "volFields.H"

#include "kineticReactions/reversibleKineticReaction/reversibleKineticReaction.H"
#include "kineticReactions/multiplicativeMonod/multiplicativeMonod.H"
#include "equilibriumReactions/waterEquilibrium/waterEquilibrium.H"
#include "equilibriumReactions/waterCarbonEquilibrium/waterCarbonEquilibrium.H"

namespace Foam
{

class reactionManager
{
private:

    fvMesh& mesh_;

    // Reference to the dictionary already read in createFields.H
    const dictionary& reactionsDict_;

    PtrList<reversibleKineticReaction> reversibleKineticReactions_;
    PtrList<multiplicativeMonod> monodReactions_;
    PtrList<waterEquilibrium> waterEquilibriumReactions_;
    PtrList<waterCarbonEquilibrium> waterCarbonEquilibriumReactions_;

    void readKineticReactions();
    void readEquilibriumReactions();

    tmp<volScalarField> zeroSource
    (
        const word& speciesName
    ) const;

public:

    reactionManager
    (
        fvMesh& mesh,
        const dictionary& reactionsDict
    );

    tmp<volScalarField> kineticSource
    (
        const word& speciesName
    ) const;

    void applyEquilibriumReactions();

    label nKineticReactions() const
    {
        return monodReactions_.size() 
        + reversibleKineticReactions_.size()
        + waterCarbonEquilibriumReactions_.size();
    }

    label nEquilibriumReactions() const
    {
        return
            waterEquilibriumReactions_.size()
          + waterCarbonEquilibriumReactions_.size();
    }

    // Accessors for waterEquilibrium diagnostics  
    bool hasWaterEquilibrium() const
    {
        return !waterEquilibriumReactions_.empty();
    }

    label nWaterEquilibriumReactions() const
    {
        return waterEquilibriumReactions_.size();
    }

    const waterEquilibrium& waterEquilibriumReaction
    (
        const label reactioni
    ) const
    {
        return waterEquilibriumReactions_[reactioni];
    }

    bool hasWaterCarbonEquilibrium() const
    {
        return !waterCarbonEquilibriumReactions_.empty();
    }
};

} // End namespace Foam

#endif

