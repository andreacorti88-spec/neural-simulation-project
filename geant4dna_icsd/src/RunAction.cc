#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

RunAction::RunAction() : G4UserRunAction() {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  // Filename encodes the run number; rename/move the file after each
  // run if you script over several energies/options (see
  // macros/scan_energies.mac and README.md for the recommended
  // one-run-per-configuration workflow).
  std::ostringstream fname;
  fname << "icsd_run" << run->GetRunID() << ".csv";
  fOutFile.open(fname.str());
  fOutFile << "# ionizations_in_target_per_event\n";
}

void RunAction::EndOfRunAction(const G4Run*)
{
  if (fOutFile.is_open()) fOutFile.close();
}

void RunAction::RecordEvent(G4int ionizationCount)
{
  if (fOutFile.is_open()) fOutFile << ionizationCount << "\n";
}
