#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(G4int particleType, G4double zStart)
  : G4VUserPrimaryGeneratorAction(),
    fParticleType(particleType),
    fSampleRadius(5.*CLHEP::mm),     // v2: matches a real CNAO clinical spot
                                        // sigma (several mm) -- see README,
                                        // sezione sulla scoperta della
                                        // diffusione laterale multipla
    fZStart(zStart)
{
  fParticleGun = new G4ParticleGun(1);

  if (fParticleType != 1) {
    // Proton: G4ParticleTable is ready to use immediately, unlike
    // G4IonTable (see GeneratePrimaries() for why carbon-12 is set
    // up lazily instead of here).
    auto particleTable = G4ParticleTable::GetParticleTable();
    auto proton = particleTable->FindParticle("proton");
    fParticleGun->SetParticleDefinition(proton);
    fParticleGun->SetParticleEnergy(150.*CLHEP::MeV);  // overridden by /gun/energy
  }
  // For carbon-12 (particleType==1), the particle definition is left
  // unset here on purpose.

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Carbon-12 must be looked up via G4IonTable, which is only fully
  // populated (GenericIon ready) once the run has actually started --
  // doing this in the constructor throws "GenericIon is not ready".
  // Done once, lazily, on the first call instead.
  if (fParticleType == 1 && !fIonReady) {
    // NOTE: deliberately does NOT call SetParticleEnergy here -- main()
    // already issued /gun/energy (via the UI, before /run/beamOn) and
    // that value must survive; SetParticleDefinition/SetParticleCharge
    // do not touch the stored energy.
    G4ParticleDefinition* carbon12 =
        G4IonTable::GetIonTable()->GetIon(6, 12, 0.*CLHEP::keV);
    if (!carbon12) {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "Ion001",
                   FatalException, "G4IonTable::GetIon(6,12,0) returned null");
    }
    fParticleGun->SetParticleDefinition(carbon12);
    fParticleGun->SetParticleCharge(6.*CLHEP::eplus);
    fIonReady = true;
  }

  // Uniform sampling on a disk: r = R*sqrt(u), theta = 2*pi*v
  G4double r = fSampleRadius * std::sqrt(G4UniformRand());
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4double x = r * std::cos(phi);
  G4double y = r * std::sin(phi);

  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, fZStart));
  fParticleGun->GeneratePrimaryVertex(event);
}
