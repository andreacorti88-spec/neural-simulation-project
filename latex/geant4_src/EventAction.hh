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

    void AddEdep(G4int neuronIndex, G4double edep, G4double stepLength);
    void MarkPrimaryHit(G4int neuronIndex);

  private:
    RunAction* fRunAction = nullptr;
    std::array<G4double, N_NEURONS> fEdep{};
    // Numeratore del LET dose-mediato, Sum(dE_i^2/dx_i) per neurone
    // (Kanai et al. -- lo stesso concetto usato a NIRS/HIMAC per il
    // modello RBE clinico): pesare ogni step per il proprio deposito
    // di energia, non contarlo una volta sola come nel LET track-medio.
    std::array<G4double, N_NEURONS> fEdepSqOverDx{};
    std::array<G4bool, N_NEURONS> fPrimaryHit{};
};

#endif
