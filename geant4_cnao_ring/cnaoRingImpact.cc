// cnaoRingImpact -- physical dose/hit map of a CNAO-style clinical
// ion beam (proton or fully-stripped carbon-12) across a ring of
// 100 tissue-equivalent "neuron" spheres, the physical counterpart
// of the N_E=100 excitatory ring in the spiking ring attractor of
// neuroni-progetto (dialogo_spiking_*.py, sezioni 5.6-5.20).
//
// Physics: QBBC (Geant4's reference physics list for clinical
// proton/ion energies -- NOT Geant4-DNA, which is only valid up to
// ~1 MeV and has no validated ion transport model at these
// energies; see PhysicsList.hh and README.md).
//
// Usage:
//   ./cnaoRingImpact [macro] [particleType] [energy_MeV] [ringDepth_mm]
//
//   macro          path to a .mac file (default: macros/run.mac)
//   particleType   0 = proton, 1 = carbon-12 (default 0)
//   energy_MeV     KINETIC energy of the primary, TOTAL (not per
//                  nucleon) for carbon-12 (default: 150 for proton,
//                  2400 for carbon-12, i.e. 200 MeV/u -- both mid-
//                  range clinical CNAO values)
//   ringDepth_mm   distance from the beam entry face to the ring
//                  plane (default: 2mm, the shallow/plateau-only
//                  geometry of sections 5.21-5.24). Use braggProfile/
//                  to find the Bragg peak depth for a given energy
//                  and pass it here to test the ring AT the peak --
//                  the clinically relevant scenario.
//
// Example, 150 MeV protons, ring in the shallow plateau (original):
//   ./cnaoRingImpact macros/run.mac 0 150 2
// Example, 70 MeV protons, ring AT the Bragg peak (39.5mm, found with
// braggProfile):
//   ./cnaoRingImpact macros/run.mac 0 70 39.5

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <ctime>
#include <string>

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv)
{
  G4long seed = static_cast<G4long>(std::time(nullptr));
  CLHEP::HepRandom::setTheSeed(seed);

  G4String macro = (argc > 1) ? argv[1] : "macros/run.mac";
  G4int particleType = (argc > 2) ? std::stoi(argv[2]) : 0;
  G4double energyMeV = (argc > 3) ? std::stod(argv[3])
                                   : (particleType == 1 ? 2400.0 : 150.0);
  G4double ringDepthMM = (argc > 4) ? std::stod(argv[4]) : 2.0;

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  auto* detector = new DetectorConstruction(ringDepthMM * CLHEP::mm);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization(particleType, detector->GetZStart()));

  runManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/gun/energy " + std::to_string(energyMeV) + " MeV");
  UImanager->ApplyCommand("/control/execute " + macro);

  delete runManager;
  return 0;
}
