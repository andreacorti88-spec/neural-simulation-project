// braggProfile -- caratterizza la curva di Bragg di un fascio CNAO
// (protone o carbonio-12) in acqua, per trovare a quale profondita'
// cade il picco a una data energia. Serve a posizionare correttamente
// l'anello di neuroni di cnaoRingImpact AL PICCO invece che nella sola
// regione di plateau (unica testata nelle sezioni 5.21-5.24 del
// resoconto), lo scenario clinicamente rilevante: al CNAO l'energia si
// sceglie apposta perche' il picco cada sul bersaglio.
//
// Geometria: fantoccio d'acqua diviso in fette sottili lungo Z (default
// 1mm x 350 = 35cm, copre l'intero range clinico dei protoni CNAO fino
// a 226.7 MeV). Fascio: pencil beam on-axis, nessuna dispersione
// laterale (basta un profilo 1D di dose in profondita').
//
// Uso: ./braggProfile [particleType 0=p,1=C12] [energy_MeV] [nEvents]
// Esempio, protone 150 MeV, 3000 primari:
//   ./braggProfile 0 150 3000

#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4VUserActionInitialization.hh"
#include "QBBC.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>

static const G4int N_SLABS = 350;
static const G4double SLAB_THICK = 1.0 * CLHEP::mm;

// ---------------- Detector ----------------
class Detector : public G4VUserDetectorConstruction {
  public:
    G4VPhysicalVolume* Construct() override {
      auto nist = G4NistManager::Instance();
      auto water = nist->FindOrBuildMaterial("G4_WATER");
      G4double totalZ = N_SLABS * SLAB_THICK;
      auto solidWorld = new G4Box("World", 5*CLHEP::cm, 5*CLHEP::cm, 0.6*totalZ);
      auto logicWorld = new G4LogicalVolume(solidWorld, water, "World");
      auto physWorld = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0, true);

      auto solidSlab = new G4Box("Slab", 5*CLHEP::cm, 5*CLHEP::cm, 0.5*SLAB_THICK);
      auto logicSlab = new G4LogicalVolume(solidSlab, water, "Slab");
      for (G4int i = 0; i < N_SLABS; ++i) {
        G4double z = -0.5*totalZ + (i + 0.5)*SLAB_THICK;
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,z), logicSlab, "Slab",
                           logicWorld, false, i, true);
      }
      return physWorld;
    }
};

// ---------------- Primary generator ----------------
class Primary : public G4VUserPrimaryGeneratorAction {
  public:
    Primary(G4int type) : fType(type) {
      fGun = new G4ParticleGun(1);
      if (fType != 1) {
        auto proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
        fGun->SetParticleDefinition(proton);
      }
      fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
      G4double totalZ = N_SLABS * SLAB_THICK;
      fZStart = -0.5*totalZ - 1*CLHEP::mm;
    }
    ~Primary() override { delete fGun; }

    void GeneratePrimaries(G4Event* evt) override {
      if (fType == 1 && !fIonReady) {
        auto c12 = G4IonTable::GetIonTable()->GetIon(6, 12, 0.*CLHEP::keV);
        fGun->SetParticleDefinition(c12);
        fGun->SetParticleCharge(6.*CLHEP::eplus);
        fIonReady = true;
      }
      fGun->SetParticlePosition(G4ThreeVector(0,0,fZStart));
      fGun->GeneratePrimaryVertex(evt);
    }

  private:
    G4ParticleGun* fGun;
    G4int fType;
    G4bool fIonReady = false;
    G4double fZStart;
};

// ---------------- Run action (accumulates + writes CSV) ----------------
class RunAct : public G4UserRunAction {
  public:
    std::vector<G4double> sumEdep;
    RunAct() : sumEdep(N_SLABS, 0.) {}

    void BeginOfRunAction(const G4Run*) override { std::fill(sumEdep.begin(), sumEdep.end(), 0.); }

    void AddEdep(G4int slab, G4double edep) {
      if (slab >= 0 && slab < N_SLABS) sumEdep[slab] += edep;
    }

    void EndOfRunAction(const G4Run* run) override {
      std::ostringstream fname;
      fname << "bragg_run" << run->GetRunID() << ".csv";
      std::ofstream out(fname.str());
      out << "slab_index,depth_mm,sum_edep_MeV\n";
      for (G4int i = 0; i < N_SLABS; ++i) {
        G4double depth = (i + 0.5) * SLAB_THICK / CLHEP::mm;
        out << i << "," << depth << "," << sumEdep[i]/CLHEP::MeV << "\n";
      }
    }
};

// ---------------- Stepping action ----------------
class SteppingAct : public G4UserSteppingAction {
  public:
    explicit SteppingAct(RunAct* run) : fRun(run) {}
    void UserSteppingAction(const G4Step* step) override {
      auto vol = step->GetPreStepPoint()->GetTouchable()->GetVolume();
      if (!vol || vol->GetName() != "Slab") return;
      G4int slab = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
      G4double edep = step->GetTotalEnergyDeposit();
      if (edep > 0.) fRun->AddEdep(slab, edep);
    }
  private:
    RunAct* fRun;
};

// ---------------- Action initialization ----------------
class ActionInit : public G4VUserActionInitialization {
  public:
    explicit ActionInit(G4int type) : fType(type) {}
    void BuildForMaster() const override { SetUserAction(new RunAct()); }
    void Build() const override {
      SetUserAction(new Primary(fType));
      auto run = new RunAct();
      SetUserAction(run);
      SetUserAction(new SteppingAct(run));
    }
  private:
    G4int fType;
};

int main(int argc, char** argv) {
  CLHEP::HepRandom::setTheSeed(static_cast<G4long>(std::time(nullptr)));

  G4int type = (argc > 1) ? std::stoi(argv[1]) : 0;
  G4double energyMeV = (argc > 2) ? std::stod(argv[2]) : 150.0;
  G4int nEvents = (argc > 3) ? std::stoi(argv[3]) : 3000;

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new Detector());
  runManager->SetUserInitialization(new QBBC(1));
  runManager->SetUserInitialization(new ActionInit(type));
  runManager->Initialize();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/gun/energy " + std::to_string(energyMeV) + " MeV");
  ui->ApplyCommand("/run/beamOn " + std::to_string(nEvents));

  delete runManager;
  return 0;
}
