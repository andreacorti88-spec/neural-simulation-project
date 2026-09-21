#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <array>
#include <vector>

class RunAction;

// Per-event accumulator: for each of the 100 ring positions, the
// total energy deposited (all tracks, primary + secondaries) and
// whether the PRIMARY ion/proton track itself crossed that neuron's
// volume at least once this event. Flushed into RunAction's running
// totals at the end of every event.
class EventAction : public G4UserEventAction
{
  public:
    static const G4int N_NEURONS = 100;

    explicit EventAction(RunAction* runAction);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddEdep(G4int neuronIndex, G4double edep, G4double stepLength);
    void MarkPrimaryHit(G4int neuronIndex);
    void AddSecondaryElectron(G4int neuronIndex, G4double kineticEnergy);

  private:
    RunAction* fRunAction = nullptr;
    std::array<G4double, N_NEURONS> fEdep{};
    // Numeratore del LET dose-mediato, Sum(dE_i^2/dx_i) per neurone
    // (Kanai et al. -- lo stesso concetto usato a NIRS/HIMAC per il
    // modello RBE clinico): pesare ogni step per il proprio deposito
    // di energia, non contarlo una volta sola come nel LET track-medio.
    std::array<G4double, N_NEURONS> fEdepSqOverDx{};
    std::array<G4bool, N_NEURONS> fPrimaryHit{};
    // Energia cinetica (MeV) di ogni elettrone secondario NATO dentro
    // un neurone in questo evento -- lo spettro necessario per valutare
    // se un accoppiamento con Geant4-DNA (valido fino a ~1 MeV) ha
    // senso fisico, sezione 5.39.
    std::vector<G4double> fSecondaryElectronEnergies{};
    // Conteggio di elettroni secondari PER NEURONE (sezione 5.44 --
    // finora la 5.39 aggregava lo spettro su tutto l'anello, senza
    // sapere da quale neurone venisse ciascun elettrone).
    std::array<G4long, N_NEURONS> fSecondaryElectronCount{};
    // Danno DSB atteso PER NEURONE, stimato assegnando a ogni elettrone
    // secondario il tasso DSB/primario del punto rappresentativo piu'
    // vicino tra i 6 della sezione 5.43 (approssimazione dichiarata:
    // non e' una nuova simulazione di danno per neurone, riusa i tassi
    // gia' validati).
    std::array<G4double, N_NEURONS> fEstimatedDSB{};
};

#endif
