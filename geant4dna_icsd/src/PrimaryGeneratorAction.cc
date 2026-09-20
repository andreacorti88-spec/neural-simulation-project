#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun(1);
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto electron = particleTable->FindParticle("e-");
  fParticleGun->SetParticleDefinition(electron);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleGun->SetParticleEnergy(100.*CLHEP::eV);   // overridden by /gun/energy
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Isotropic direction: sample uniformly on the unit sphere.
  G4double cosTheta = 2.*G4UniformRand() - 1.;
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4ThreeVector dir(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleGun->SetParticleMomentumDirection(dir);
  fParticleGun->GeneratePrimaryVertex(event);
}
