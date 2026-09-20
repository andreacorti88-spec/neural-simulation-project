#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"

DetectorConstruction::DetectorConstruction(G4double ringDepthFromEntry)
  : G4VUserDetectorConstruction(),
    fRingDepth(ringDepthFromEntry),
    fRingRadius(500.*CLHEP::um),        // ~1mm diameter ring, cortical
                                          // microcircuit scale (see README)
    fSomaRadius(10.*CLHEP::um)          // ~20um soma diameter, typical
                                          // cortical pyramidal neuron (README)
{
  // World sized to comfortably contain the beam's full path from entry
  // to the ring, plus a safety margin (proton/ion range straggling and
  // the lateral 550um sampling disk).
  fWorldHalfXY = 5.*CLHEP::mm;
  fWorldHalfZ = 0.5*fRingDepth + 5.*CLHEP::mm;
  fZStart = -fWorldHalfZ + 0.5*CLHEP::mm;  // just inside the entry face
}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  // Water as tissue-equivalent surrogate, the same convention used in
  // clinical proton/ion dosimetry and in nanoICSD -- keeps this run
  // directly comparable in dose units to standard hadrontherapy
  // Monte Carlo practice.
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

  // --- World: tissue block along the beam axis, ring placed at
  //     Z = fZStart + fRingDepth (the depth requested by main(), e.g.
  //     the Bragg peak depth found with braggProfile/) ---
  auto solidWorld = new G4Box("World", fWorldHalfXY, fWorldHalfXY, fWorldHalfZ);
  auto logicWorld  = new G4LogicalVolume(solidWorld, water, "World");
  auto physWorld   = new G4PVPlacement(nullptr, {}, logicWorld, "World",
                                        nullptr, false, 0, true);

  G4double ringZ = fZStart + fRingDepth;

  // --- Ring of N_NEURONS somata in the Z=ringZ plane (the requested
  //     depth), positions matching xi = i/N_NEURONS in the Brian2
  //     model (dialogo_spiking_*.py) ---
  fLogicNeuron.reserve(N_NEURONS);
  for (G4int i = 0; i < N_NEURONS; ++i) {
    G4double xi = static_cast<G4double>(i) / N_NEURONS;
    G4double theta = CLHEP::twopi * xi;
    G4double x = fRingRadius * std::cos(theta);
    G4double y = fRingRadius * std::sin(theta);

    auto solidNeuron = new G4Orb("Neuron", fSomaRadius);
    auto logicNeuron = new G4LogicalVolume(solidNeuron, water, "Neuron");
    // Copy number IS the neuron index i -- SteppingAction reads it
    // back to know which ring position (xi) was hit.
    new G4PVPlacement(nullptr, G4ThreeVector(x, y, ringZ), logicNeuron,
                       "Neuron", logicWorld, false, i, true);
    fLogicNeuron.push_back(logicNeuron);

    auto visNeuron = new G4VisAttributes(G4Colour(0.85, 0.25, 0.2, 0.6));
    logicNeuron->SetVisAttributes(visNeuron);
  }

  auto visWorld = new G4VisAttributes(G4Colour(0.5, 0.6, 1.0, 0.03));
  visWorld->SetForceWireframe(true);
  logicWorld->SetVisAttributes(visWorld);

  // Sezione 5.39: senza un taglio di produzione dedicato, il default di
  // QBBC (~1mm) e' molto piu' grande dei neuroni (10um di raggio) --
  // quasi nessun elettrone secondario (raggio delta) ha un range
  // sufficiente per essere generato come track esplicito dentro un
  // volume cosi' piccolo; la sua energia viene invece depositata
  // localmente come perdita di energia continua, INVISIBILE a un
  // tracking per-particella (verificato: zero elettroni secondari
  // registrati su 500.000 eventi prima di questa correzione). Una
  // G4Region dedicata ai volumi "Neuron", con un taglio molto piu'
  // fine (100nm, sotto la scala del neurone), fa si' che i raggi delta
  // prodotti li' dentro vengano effettivamente tracciati come
  // particelle esplicite -- il prerequisito per qualunque accoppiamento
  // con Geant4-DNA (che opera per definizione su track espliciti).
  auto* neuronRegion = new G4Region("NeuronRegion");
  auto* neuronCuts = new G4ProductionCuts();
  neuronCuts->SetProductionCut(100.*CLHEP::nm);
  neuronRegion->SetProductionCuts(neuronCuts);
  for (auto* logicNeuron : fLogicNeuron) {
    neuronRegion->AddRootLogicalVolume(logicNeuron);
  }

  return physWorld;
}
