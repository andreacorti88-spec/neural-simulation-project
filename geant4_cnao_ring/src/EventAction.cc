#include "EventAction.hh"
#include "RunAction.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

namespace {
// Sezione 5.43: 6 punti rappresentativi dello spettro reale (percentili
// 10/25/50/75/90/99), con il tasso DSB/primario validato per ciascuno
// (moleculardna, geant4dna_ssbdsb). Riusati qui, non ricalcolati.
constexpr G4int kNBins = 6;
constexpr G4double kBinEnergyKeV[kNBins] = {1.09, 1.28, 1.86, 3.43, 7.33, 28.82};
constexpr G4double kBinDSBRate[kNBins]   = {0.0103, 0.0083, 0.0137, 0.0243, 0.0430, 0.0163};

G4double DSBRateForEnergy(G4double keV)
{
  G4int best = 0;
  G4double bestDiff = std::abs(keV - kBinEnergyKeV[0]);
  for (G4int i = 1; i < kNBins; ++i) {
    G4double diff = std::abs(keV - kBinEnergyKeV[i]);
    if (diff < bestDiff) { bestDiff = diff; best = i; }
  }
  return kBinDSBRate[best];
}
}  // namespace

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
  fSecondaryElectronCount.fill(0);
  fEstimatedDSB.fill(0.);
  fFragmentEdep.fill(0.);
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
  fSecondaryElectronCount[neuronIndex] += 1;
  fEstimatedDSB[neuronIndex] += DSBRateForEnergy(kineticEnergy / CLHEP::keV);
}

void EventAction::AddFragmentEdep(G4int neuronIndex, G4double edep)
{
  if (neuronIndex < 0 || neuronIndex >= N_NEURONS) return;
  fFragmentEdep[neuronIndex] += edep;
}

void EventAction::EndOfEventAction(const G4Event*)
{
  fRunAction->RecordEvent(fEdep, fEdepSqOverDx, fPrimaryHit, fSecondaryElectronEnergies,
                           fSecondaryElectronCount, fEstimatedDSB, fFragmentEdep);
}
