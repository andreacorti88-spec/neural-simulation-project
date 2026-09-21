// nanoICSD -- reproduces the ICSD scoring test case of
// Villagrasa, Baiocco et al., PLOS One 21(1):e0340500 (2026).
//
// Usage:
//   ./nanoICSD [macro] [physicsOption] [targetDiameter_nm] [seed]
//
//   macro              path to a .mac file (default: macros/run.mac)
//   physicsOption       2, 4 or 6 -> Geant4-DNA option 2/4/6 (default 2)
//   targetDiameter_nm    8 or 100, matching the paper's geometry (default 8)
//   seed                 explicit random seed (default: system clock, i.e.
//                        a DIFFERENT Monte Carlo realization on every run
//                        with no 4th argument -- pass an explicit seed for
//                        reproducible runs or multi-seed MC verification of
//                        M1/F2/F3, the same discipline applied to
//                        cnaoRingImpact after section 5.32 found its
//                        per-neuron dose map was seed-unstable).
//
// Example, reproducing the 100 eV / option 2 / 8 nm case:
//   ./nanoICSD macros/run.mac 2 8
// Example, same point with an explicit seed for MC verification:
//   ./nanoICSD macros/run.mac 2 8 12345
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
  G4String macro = (argc > 1) ? argv[1] : "macros/run.mac";
  G4int physicsOption = (argc > 2) ? std::stoi(argv[2]) : 2;
  G4double targetDiameter_nm = (argc > 3) ? std::stod(argv[3]) : 8.0;
  // Independent seed per process -- important if you parallelize
  // energy points / options across several batch jobs. Default: system
  // clock (different every run with no 4th argument); pass one
  // explicitly for reproducibility or multi-seed MC verification.
  G4long seed = (argc > 4) ? std::stol(argv[4])
                            : static_cast<G4long>(std::time(nullptr));
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << "Seed usato: " << seed << G4endl;

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
