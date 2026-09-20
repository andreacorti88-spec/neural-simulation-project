#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream>

class G4Run;

// Opens one plain-text file per run and writes one line per event:
// the number of ionizations recorded inside the target sphere for
// that primary electron. Post-processing (Python, see
// analysis/analyze_icsd.py) turns this file into the ICSD histogram
// and derived M1 / F2 / F3, in the same way Table 3/4 of the paper
// were built from the raw per-event ionization counts.
class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void RecordEvent(G4int ionizationCount);

  private:
    std::ofstream fOutFile;
};

#endif
