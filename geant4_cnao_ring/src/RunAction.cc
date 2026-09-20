#include "RunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

RunAction::RunAction() : G4UserRunAction() {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
  fSumEdep.fill(0.);
  fSumPrimaryHits.fill(0);
  fNEvents = 0;
}

void RunAction::RecordEvent(const std::array<G4double, N_NEURONS>& edep,
                             const std::array<G4bool, N_NEURONS>& primaryHit)
{
  for (G4int i = 0; i < N_NEURONS; ++i) {
    fSumEdep[i] += edep[i];
    if (primaryHit[i]) fSumPrimaryHits[i] += 1;
  }
  fNEvents += 1;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  std::ostringstream fname;
  fname << "neuron_dose_run" << run->GetRunID() << ".csv";
  std::ofstream out(fname.str());
  out << "# n_events=" << fNEvents << "\n";
  out << "neuron_index,xi,sum_edep_MeV,mean_edep_keV_per_event,"
         "primary_hit_count,primary_hit_fraction\n";
  for (G4int i = 0; i < N_NEURONS; ++i) {
    G4double xi = static_cast<G4double>(i) / N_NEURONS;
    G4double sumEdepMeV = fSumEdep[i] / CLHEP::MeV;
    G4double meanEdepKeV = (fNEvents > 0)
        ? (fSumEdep[i] / CLHEP::keV) / static_cast<G4double>(fNEvents) : 0.;
    G4double hitFraction = (fNEvents > 0)
        ? static_cast<G4double>(fSumPrimaryHits[i]) / static_cast<G4double>(fNEvents) : 0.;
    out << i << "," << xi << "," << sumEdepMeV << "," << meanEdepKeV << ","
        << fSumPrimaryHits[i] << "," << hitFraction << "\n";
  }
  out.close();
}
