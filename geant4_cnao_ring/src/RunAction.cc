#include "RunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

RunAction::RunAction() : G4UserRunAction() {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
  fSumEdep.fill(0.);
  fSumEdepSqOverDx.fill(0.);
  fSumPrimaryHits.fill(0);
  fNEvents = 0;
  fSecondaryElectronSpectrum.clear();
}

void RunAction::RecordEvent(const std::array<G4double, N_NEURONS>& edep,
                             const std::array<G4double, N_NEURONS>& edepSqOverDx,
                             const std::array<G4bool, N_NEURONS>& primaryHit,
                             const std::vector<G4double>& secondaryElectronEnergies)
{
  for (G4int i = 0; i < N_NEURONS; ++i) {
    fSumEdep[i] += edep[i];
    fSumEdepSqOverDx[i] += edepSqOverDx[i];
    if (primaryHit[i]) fSumPrimaryHits[i] += 1;
  }
  fNEvents += 1;
  for (G4double e : secondaryElectronEnergies) {
    fSecondaryElectronSpectrum.push_back(e);
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  std::ostringstream fname;
  fname << "neuron_dose_run" << run->GetRunID() << ".csv";
  std::ofstream out(fname.str());
  out << "# n_events=" << fNEvents << "\n";
  out << "neuron_index,xi,sum_edep_MeV,mean_edep_keV_per_event,"
         "primary_hit_count,primary_hit_fraction,dose_averaged_LET_keV_per_um\n";
  for (G4int i = 0; i < N_NEURONS; ++i) {
    G4double xi = static_cast<G4double>(i) / N_NEURONS;
    G4double sumEdepMeV = fSumEdep[i] / CLHEP::MeV;
    G4double meanEdepKeV = (fNEvents > 0)
        ? (fSumEdep[i] / CLHEP::keV) / static_cast<G4double>(fNEvents) : 0.;
    G4double hitFraction = (fNEvents > 0)
        ? static_cast<G4double>(fSumPrimaryHits[i]) / static_cast<G4double>(fNEvents) : 0.;
    // LET dose-mediato = Sum(dE_i^2/dx_i) / Sum(dE_i) (Kanai et al. 1999,
    // lo stesso concetto usato a NIRS/HIMAC per il modello RBE clinico:
    // pesa ogni step per il proprio deposito di energia).
    G4double doseAvgLET = (fSumEdep[i] > 0.)
        ? (fSumEdepSqOverDx[i] / fSumEdep[i]) / (CLHEP::keV / CLHEP::um) : 0.;
    out << i << "," << xi << "," << sumEdepMeV << "," << meanEdepKeV << ","
        << fSumPrimaryHits[i] << "," << hitFraction << "," << doseAvgLET << "\n";
  }
  out.close();

  // Spettro degli elettroni secondari nati dentro un neurone (sezione
  // 5.39): un'energia cinetica (keV) per riga, cosi' come creati --
  // input per valutare l'accoppiamento con Geant4-DNA.
  std::ostringstream specName;
  specName << "secondary_electron_spectrum_run" << run->GetRunID() << ".csv";
  std::ofstream specOut(specName.str());
  specOut << "# n_events=" << fNEvents << " n_secondary_electrons="
          << fSecondaryElectronSpectrum.size() << "\n";
  specOut << "kinetic_energy_keV\n";
  for (G4double e : fSecondaryElectronSpectrum) {
    specOut << (e / CLHEP::keV) << "\n";
  }
  specOut.close();
}
