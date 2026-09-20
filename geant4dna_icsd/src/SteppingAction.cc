#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"

SteppingAction::SteppingAction(EventAction* eventAction,
                                const G4LogicalVolume* targetVolume)
  : G4UserSteppingAction(), fEventAction(eventAction), fTargetVolume(targetVolume)
{}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4StepPoint* preStep = step->GetPreStepPoint();
  const G4VPhysicalVolume* volume = preStep->GetTouchableHandle()->GetVolume();
  if (!volume || volume->GetName() != "Target") return;

  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  if (!process) return;

  if (process->GetProcessName() == "e-_G4DNAIonisation") {
    fEventAction->AddIonization();
  }
}
