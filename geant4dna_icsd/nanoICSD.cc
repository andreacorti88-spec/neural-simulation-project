// nanoICSD -- reproduces the ICSD scoring test case of
// Villagrasa, Baiocco et al., PLOS One 21(1):e0340500 (2026).
//
// Usage:
//   ./nanoICSD [macro] [physicsOption] [targetDiameter_nm]
//
//   macro              path to a .mac file (default: macros/run.mac)
//   physicsOption       2, 4 or 6 -> Geant4-DNA option 2/4/6 (default 2)
//   targetDiameter_nm    8 or 100, matching the paper's geometry (default 8)
//
// Example, reproducing the 100 eV / option 2 / 8 nm case:
//   ./nanoICSD macros/run.mac 2 8
//
// See README.md for how to loop this over the full energy grid used
// in the paper (20, 50, 100, 300, 600, 1000, 5000, 10000 eV).

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <ctime>

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv)
{
  // Independent seed per process -- important if you parallelize
  // energy points / options across several batch jobs.
  G4long seed = static_cast<G4long>(std::time(nullptr));
  CLHEP::HepRandom::setTheSeed(seed);

  G4String macro = (argc > 1) ? argv[1] : "macros/run.mac";
  G4int physicsOption = (argc > 2) ? std::stoi(argv[2]) : 2;
  G4double targetDiameter_nm = (argc > 3) ? std::stod(argv[3]) : 8.0;

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  auto* detector = new DetectorConstruction();
  detector->SetTargetDiameter(targetDiameter_nm * CLHEP::nanometer);
  runManager->SetUserInitialization(detector);

  runManager->SetUserInitialization(new PhysicsList(physicsOption));
  runManager->SetUserInitialization(new ActionInitialization(detector));

  runManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4String command = "/control/execute ";
  UImanager->ApplyCommand(command + macro);

  delete runManager;
  return 0;
}
