#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

class G4LogicalVolume;

// Geometry: a water-equivalent tissue block (the surrounding
// neuropil/extracellular space) containing a ring of N_NEURONS
// spherical "somata" -- the physical stand-in for the N_E=100
// excitatory neurons of the ring attractor in neuroni-progetto
// (dialogo_spiking_*.py). Only the excitatory population carries a
// meaningful ring position x in that model, so only it is given a
// physical counterpart here.
//
// Ring radius and soma radius are literature-typical values for a
// cortical microcircuit (see README.md for the specific sources),
// not values read off any specific measurement of the simulated
// network -- the ring attractor itself is a topological, not a
// metric, model.
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    static const G4int N_NEURONS = 100;

    // ringDepthFromEntry: distance along Z from the beam entry face to
    // the ring plane (e.g. the Bragg peak depth for a given beam --
    // see braggProfile/). Default (2mm) matches the original shallow/
    // plateau-only geometry of sections 5.21-5.24.
    explicit DetectorConstruction(G4double ringDepthFromEntry = 2.*CLHEP::mm);
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    const std::vector<G4LogicalVolume*>& GetNeuronVolumes() const { return fLogicNeuron; }
    G4double GetRingRadius() const { return fRingRadius; }
    G4double GetSomaRadius() const { return fSomaRadius; }
    G4double GetZStart() const { return fZStart; }

  private:
    G4double fRingDepth;
    G4double fWorldHalfZ;
    G4double fWorldHalfXY;
    G4double fZStart;
    G4double fRingRadius;
    G4double fSomaRadius;
    std::vector<G4LogicalVolume*> fLogicNeuron;
};

#endif
