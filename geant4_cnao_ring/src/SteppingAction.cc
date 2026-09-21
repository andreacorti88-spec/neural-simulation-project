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
#include "G4RunManager.hh"
#include "G4Event.hh"
#include <cstdlib>

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
  if (edep > 0.) fEventAction->AddEdep(neuronIndex, edep, step->GetStepLength());

  // Diagnostica temporanea (sezione 5.36): stampa ogni step con
  // deposito >0 nel neurone indicato da CNAO_DEBUG_NEURON, per
  // identificare la particella/processo responsabile di un outlier
  // di LET -- stesso pattern usato per il bug del geantino.
  static const char* dbgEnv = std::getenv("CNAO_DEBUG_NEURON");
  if (dbgEnv && neuronIndex == std::atoi(dbgEnv) && edep > 0.) {
    const G4Track* track = step->GetTrack();
    const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4cout << "[CNAO_DEBUG_NEURON] event=" << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
           << " particella=" << track->GetParticleDefinition()->GetParticleName()
           << " parentID=" << track->GetParentID()
           << " trackID=" << track->GetTrackID()
           << " processo=" << (proc ? proc->GetProcessName() : "N/A")
           << " edep_keV=" << edep/CLHEP::keV
           << " stepLength_um=" << step->GetStepLength()/CLHEP::um
           << " KE_MeV=" << track->GetKineticEnergy()/CLHEP::MeV
           << G4endl;
  }

  // ParentID==0 identifies the primary proton/ion track itself
  // (as opposed to any secondary electron, delta-ray, or nuclear
  // fragment produced along the way).
  const G4Track* track = step->GetTrack();
  if (track->GetParentID() == 0) {
    fEventAction->MarkPrimaryHit(neuronIndex);
  }

  // Sezione 5.39: registra l'energia cinetica di ogni elettrone
  // secondario NATO dentro questo neurone (primo step del track,
  // pre-step point ancora nel volume "Neuron") -- lo spettro
  // necessario per valutare un accoppiamento con Geant4-DNA (valido
  // solo fino a ~1 MeV), invece di assumerlo senza verificarlo.
  if (track->GetParentID() != 0
      && track->GetParticleDefinition()->GetParticleName() == "e-"
      && track->GetCurrentStepNumber() == 1) {
    fEventAction->AddSecondaryElectron(neuronIndex, preStep->GetKineticEnergy());
  }

  // Sezione 5.46: dose da FRAMMENTI NUCLEARI secondari carichi --
  // qualunque traccia non primaria (parentID!=0), non un elettrone,
  // con carica positiva (protoni, alfa, deutoni, tritoni, o ioni piu'
  // pesanti prodotti da reazioni inelastiche di un primario di
  // carbonio-12). La coda di frammentazione oltre il picco di Bragg
  // del carbonio e' una differenza fisica reale e clinicamente nota
  // tra terapia a ioni di carbonio e protonterapia (letteratura
  // Chiba/NIRS-HIMAC): questo accumulatore la rende misurabile.
  if (edep > 0. && track->GetParentID() != 0
      && track->GetParticleDefinition()->GetParticleName() != "e-"
      && track->GetParticleDefinition()->GetPDGCharge() > 0.) {
    fEventAction->AddFragmentEdep(neuronIndex, edep);
  }
}
