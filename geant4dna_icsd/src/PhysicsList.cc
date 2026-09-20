#include "PhysicsList.hh"

#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList(G4int option) : G4VModularPhysicsList()
{
  SetVerboseLevel(1);

  // Geant4-DNA physics is only valid down to a few eV and only for
  // its supported materials (liquid water here) -- this mirrors the
  // configuration used in every standard Geant4-DNA example.
  G4EmParameters::Instance()->SetApplyCuts(true);
  G4EmParameters::Instance()->SetMinEnergy(10*CLHEP::eV);
  G4EmParameters::Instance()->SetMaxEnergy(1*CLHEP::MeV);
  G4EmParameters::Instance()->SetBuildCSDARange(true);

  switch (option) {
    case 2:
      RegisterPhysics(new G4EmDNAPhysics_option2());
      G4cout << "PhysicsList: using G4EmDNAPhysics_option2 (Geant4-DNA opt.2)" << G4endl;
      break;
    case 4:
      RegisterPhysics(new G4EmDNAPhysics_option4());
      G4cout << "PhysicsList: using G4EmDNAPhysics_option4 (Geant4-DNA opt.4)" << G4endl;
      break;
    case 6:
      RegisterPhysics(new G4EmDNAPhysics_option6());
      G4cout << "PhysicsList: using G4EmDNAPhysics_option6 (Geant4-DNA opt.6)" << G4endl;
      break;
    default:
      RegisterPhysics(new G4EmDNAPhysics());
      G4cout << "PhysicsList: using default G4EmDNAPhysics" << G4endl;
      break;
  }
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{
  // Geant4-DNA is a discrete, cut-free physics regime by design;
  // this call just keeps SetCuts() well-defined for the base class.
  SetCutsWithDefault();
}
