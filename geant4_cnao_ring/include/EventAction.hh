#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <array>

class RunAction;

// Per-event accumulator: for each of the 100 ring positions, the
// total energy deposited (all tracks, primary + secondaries) and
// whether the PRIMARY ion/proton track itself crossed that neuron's
// volume at least once this event. Flushed into RunAction's running
// totals at the end of every event.
class EventAction : public G4UserEventAction
{
  public:
    static const G4int N_NEURONS = 100;

    explicit EventAction(RunAction* runAction);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddEdep(G4int neuronIndex, G4double edep);
    void MarkPrimaryHit(G4int neuronIndex);

  private:
    RunAction* fRunAction = nullptr;
    std::array<G4double, N_NEURONS> fEdep{};
    std::array<G4bool, N_NEURONS> fPrimaryHit{};
};

#endif
