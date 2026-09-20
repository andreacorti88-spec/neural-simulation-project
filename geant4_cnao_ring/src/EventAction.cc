#include "EventAction.hh"
#include "RunAction.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(), fRunAction(runAction)
{}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep.fill(0.);
  fPrimaryHit.fill(false);
}

void EventAction::AddEdep(G4int neuronIndex, G4double edep)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fEdep[neuronIndex] += edep;
}

void EventAction::MarkPrimaryHit(G4int neuronIndex)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fPrimaryHit[neuronIndex] = true;
}

void EventAction::EndOfEventAction(const G4Event*)
{
  fRunAction->RecordEvent(fEdep, fPrimaryHit);
}
