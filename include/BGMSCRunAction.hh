#ifndef BGMSCRUNACTION_H
#define BGMSCRUNACTION_H

#include "G4UserRunAction.hh"
#include "BGMSCDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class BGMSCRunAction : public G4UserRunAction
{
public:

    BGMSCRunAction();
    virtual ~BGMSCRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
};

#endif


