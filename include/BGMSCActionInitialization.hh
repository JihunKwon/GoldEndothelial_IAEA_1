#ifndef BGMSCACTIONINITIASLIZATION_HH
#define BGMSCACTIONINITIASLIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "BGMSCDetectorConstruction.hh"

class BGMSCActionInitialization : public G4VUserActionInitialization
{
  public:

    BGMSCActionInitialization(BGMSCDetectorConstruction* geometry);
    virtual ~BGMSCActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif // BGMSCACTIONINITIASLIZATION_HH
