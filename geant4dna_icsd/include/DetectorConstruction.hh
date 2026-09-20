#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;

// Geometry: a large liquid-water "world" sphere containing a small
// liquid-water target sphere at its centre (the nanometric scoring
// volume). The target diameter matches the paper: 8 nm for electron
// energies up to 1 keV, 100 nm for 5 and 10 keV (set via macro).
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    // Called from the messenger / macro before Construct(), e.g.
    // "/nanoICSD/targetDiameter 8 nm"
    void SetTargetDiameter(G4double d) { fTargetDiameter = d; }
    G4double GetTargetDiameter() const { return fTargetDiameter; }

    G4LogicalVolume* GetTargetLogicalVolume() const { return fLogicTarget; }

  private:
    G4double fTargetDiameter;   // default set in .cc, override via macro
    G4double fWorldDiameter;    // kept well above the electron range studied
    G4LogicalVolume* fLogicTarget = nullptr;
};

#endif
