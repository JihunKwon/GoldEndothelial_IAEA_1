#include "BGMSCPrimaryGeneratorAction.hh"
#include "BGMSCRunAction.hh"
#include "BGMSCEventAction.hh"
#include "BGMSCSteppingAction.hh"
#include "BGMSCActionInitialization.hh"

BGMSCActionInitialization::BGMSCActionInitialization(BGMSCDetectorConstruction *geometry)
    : G4VUserActionInitialization()
{}

BGMSCActionInitialization::~BGMSCActionInitialization()
{}

void BGMSCActionInitialization::BuildForMaster() const
{
    SetUserAction(new BGMSCRunAction);
}

void BGMSCActionInitialization::Build() const
{
    SetUserAction(new BGMSCPrimaryGeneratorAction);
    SetUserAction(new BGMSCEventAction);
    SetUserAction(new BGMSCSteppingAction);
    SetUserAction(new BGMSCRunAction);
}
