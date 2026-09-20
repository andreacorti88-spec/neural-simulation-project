#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4LogicalVolume;

// Flags every step whose defining process is a Geant4-DNA ionization
// process AND whose step occurred inside the target volume, and adds
// one count to the running per-event ICSD tally. This mirrors the
// physical definition used in the paper: an "ionization" contributing
// to the ICSD is any ionizing interaction (of the primary electron
// or of any of its secondaries) that takes place within the
// nanometric target sphere.
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction, const G4LogicalVolume* targetVolume);
    ~SteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;

  private:
    EventAction* fEventAction = nullptr;
    const G4LogicalVolume* fTargetVolume = nullptr;
};

#endif
