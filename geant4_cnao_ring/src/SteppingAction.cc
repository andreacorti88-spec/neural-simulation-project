#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHandle.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
  : G4UserSteppingAction(), fEventAction(eventAction)
{}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4StepPoint* preStep = step->GetPreStepPoint();
  const G4VTouchable* touchable = preStep->GetTouchable();
  const G4VPhysicalVolume* volume = touchable->GetVolume();
  if (!volume || volume->GetName() != "Neuron") return;

  const G4int neuronIndex = touchable->GetCopyNumber();

  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) fEventAction->AddEdep(neuronIndex, edep);

  // ParentID==0 identifies the primary proton/ion track itself
  // (as opposed to any secondary electron, delta-ray, or nuclear
  // fragment produced along the way).
  if (step->GetTrack()->GetParentID() == 0) {
    fEventAction->MarkPrimaryHit(neuronIndex);
  }
}
