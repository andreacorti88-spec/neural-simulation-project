#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;

// Fires a single electron from the origin (centre of the target
// sphere) with isotropic direction, matching "a point source of
// mono-energetic electrons in its centre" (paper, Methods).
// Energy is set from the macro command /gun/energy, as usual.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;
    G4ParticleGun* GetParticleGun() { return fParticleGun; }

  private:
    G4ParticleGun* fParticleGun = nullptr;
};

#endif
