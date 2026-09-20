#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

// Pencil-beam primary: proton or fully-stripped carbon-12 ion (the
// two CNAO clinical modalities), entering along +Z from outside the
// world volume. Lateral (x,y) start position is sampled uniformly
// inside a disk that covers the whole neuron ring plus a margin --
// this approximates, at the microscopic scale simulated here, the
// locally flat portion of a real CNAO spot-scanning pencil beam
// (clinical spot sigma is several mm, i.e. much wider than the
// ~1mm ring), so the beam is effectively uniform across the ring
// on this scale; any neuron-to-neuron variation seen in the results
// is genuine microdosimetric/track-structure stochasticity, not an
// artefact of beam shape. See README.md.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    // particleType: 0 = proton, 1 = carbon-12 (fully stripped, q=+6e)
    // zStart: entry-face Z coordinate, matched to DetectorConstruction's
    // world size (which scales with the requested ring depth)
    PrimaryGeneratorAction(G4int particleType, G4double zStart);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

  private:
    G4ParticleGun* fParticleGun = nullptr;
    G4int fParticleType;
    G4bool fIonReady = false;   // G4ParticleGun defaults to "geantino"
                                  // (never nullptr), so GetParticleDefinition()
                                  // cannot be used to detect "not yet set"
    G4double fSampleRadius;
    G4double fZStart;
};

#endif
