#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

// Thin wrapper that registers one of the Geant4-DNA physics
// constructors. The "option" argument maps to the same naming used
// in the paper's Table 1 (Geant4-DNA options 2, 4, 6), so results
// can be compared line-by-line against Table 3/4/6/7 of the article.
//
//   option = 2 -> G4EmDNAPhysics_option2  (Champion elastic model)
//   option = 4 -> G4EmDNAPhysics_option4
//   option = 6 -> G4EmDNAPhysics_option6  (CPA100-based, relativistic BEB)
//   anything else -> G4EmDNAPhysics()     (Geant4-DNA default constructor)
//
// NOTE: this only reproduces the *native* cross sections shipped
// with each option. Injecting the paper's "common cross section"
// dataset instead requires supplying custom G4VEmModel / cross
// section tables -- see README.md, section "Sezioni d'urto comuni",
// for what that involves and why it is out of scope for a first
// build-and-check pass.
class PhysicsList : public G4VModularPhysicsList
{
  public:
    explicit PhysicsList(G4int option = 2);
    ~PhysicsList() override;

    void SetCuts() override;
};

#endif
