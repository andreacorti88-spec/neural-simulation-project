#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "QBBC.hh"

// QBBC is Geant4's general-purpose reference physics list, explicitly
// recommended (and used in the official "Hadrontherapy" example) for
// proton/ion beams at clinical energies (tens to hundreds of MeV) --
// a completely different regime from the Geant4-DNA track-structure
// physics used in nanoICSD, which is only valid up to ~1 MeV and only
// for electrons (see PhysicsList.cc in that project). Geant4-DNA
// track-structure physics has no validated ion transport model at
// CNAO clinical energies, so it cannot be reused here for the primary
// beam -- see README.md for the full reasoning.
class PhysicsList : public QBBC
{
  public:
    PhysicsList();
    ~PhysicsList() override;
};

#endif
