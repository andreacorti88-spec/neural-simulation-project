#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

// Counts ionizations occurring inside the target sphere during one
// primary electron's full history (primary + all secondaries down
// to the tracking cut). At the end of the event, this single
// integer is one sample of the Ionization Cluster Size Distribution
// (ICSD), exactly as scored in the paper.
class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(RunAction* runAction);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddIonization() { fIonizationCount++; }

  private:
    RunAction* fRunAction = nullptr;
    G4int fIonizationCount = 0;
};

#endif
