#include "EventAction.hh"
#include "RunAction.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(), fRunAction(runAction)
{}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep.fill(0.);
  fEdepSqOverDx.fill(0.);
  fPrimaryHit.fill(false);
  fSecondaryElectronEnergies.clear();
}

void EventAction::AddEdep(G4int neuronIndex, G4double edep, G4double stepLength)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fEdep[neuronIndex] += edep;
  if (stepLength > 0.) fEdepSqOverDx[neuronIndex] += edep * edep / stepLength;
}

void EventAction::MarkPrimaryHit(G4int neuronIndex)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fPrimaryHit[neuronIndex] = true;
}

void EventAction::AddSecondaryElectron(G4int neuronIndex, G4double kineticEnergy)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fSecondaryElectronEnergies.push_back(kineticEnergy);
}

void EventAction::EndOfEventAction(const G4Event*)
{
  fRunAction->RecordEvent(fEdep, fEdepSqOverDx, fPrimaryHit, fSecondaryElectronEnergies);
}
