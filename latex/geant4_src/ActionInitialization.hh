#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(G4int particleType, G4double zStart);
    ~ActionInitialization() override;

    void BuildForMaster() const override;
    void Build() const override;

  private:
    G4int fParticleType;
    G4double fZStart;
};

#endif
