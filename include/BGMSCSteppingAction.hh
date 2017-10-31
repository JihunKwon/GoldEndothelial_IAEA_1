#ifndef BGMSCSTEPPINGACTION_H
#define BGMSCSTEPPINGACTION_H

#include "G4UserSteppingAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

class PrimaryGeneratorAction;

class BGMSCSteppingAction : public G4UserSteppingAction {

public:
    BGMSCSteppingAction();
    ~BGMSCSteppingAction();

    void UserSteppingAction(const G4Step* aStep);

private:
    G4double EdepEndotherialTotal;
    G4double EdepEndotherial[100];
};

#endif
