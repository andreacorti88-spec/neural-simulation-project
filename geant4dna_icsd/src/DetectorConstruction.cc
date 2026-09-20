#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fTargetDiameter(8.*CLHEP::nanometer),   // paper default for <=1 keV electrons
    fWorldDiameter(2.*CLHEP::um)            // comfortably larger than any electron range studied
{}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  // G4_WATER is required (and sufficient) for the Geant4-DNA physics
  // constructors used in PhysicsList -- do not substitute a custom
  // material here, the DNA models are validated against this NIST entry.
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

  // --- World: liquid water sphere, large enough that no scored
  //     electron leaves it before slowing below the transport cut ---
  auto solidWorld = new G4Orb("World", 0.5*fWorldDiameter);
  auto logicWorld  = new G4LogicalVolume(solidWorld, water, "World");
  auto physWorld   = new G4PVPlacement(nullptr, {}, logicWorld, "World",
                                        nullptr, false, 0, true);

  // --- Target: the nanometric scoring sphere at the centre, where
  //     the primary electron starts. Ionizations occurring inside
  //     this volume are what SteppingAction accumulates into the
  //     per-event ICSD count. ---
  auto solidTarget = new G4Orb("Target", 0.5*fTargetDiameter);
  fLogicTarget = new G4LogicalVolume(solidTarget, water, "Target");
  new G4PVPlacement(nullptr, {}, fLogicTarget, "Target",
                     logicWorld, false, 0, true);

  // Visualisation only -- irrelevant for batch runs
  auto visWorld = new G4VisAttributes(G4Colour(1,1,1,0.05));
  visWorld->SetForceWireframe(true);
  logicWorld->SetVisAttributes(visWorld);
  auto visTarget = new G4VisAttributes(G4Colour(0.12,0.31,0.37,0.4));
  fLogicTarget->SetVisAttributes(visTarget);

  return physWorld;
}
