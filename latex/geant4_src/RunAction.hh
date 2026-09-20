#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <array>
#include <fstream>
#include <vector>

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    static const G4int N_NEURONS = 100;

    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void RecordEvent(const std::array<G4double, N_NEURONS>& edep,
                      const std::array<G4double, N_NEURONS>& edepSqOverDx,
                      const std::array<G4bool, N_NEURONS>& primaryHit,
                      const std::vector<G4double>& secondaryElectronEnergies);

  private:
    std::array<G4double, N_NEURONS> fSumEdep{};
    std::array<G4double, N_NEURONS> fSumEdepSqOverDx{};
    std::array<G4long, N_NEURONS> fSumPrimaryHits{};
    G4long fNEvents = 0;
    std::ofstream fOutFile;
    // Spettro di tutti gli elettroni secondari nati dentro un neurone,
    // accumulato su tutto il run (sezione 5.39).
    std::vector<G4double> fSecondaryElectronSpectrum{};
};

#endif
