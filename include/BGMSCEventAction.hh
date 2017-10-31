#ifndef BGMSCEVENTACTION_H
#define BGMSCEVENTACTION_H

#include "G4UserEventAction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class BGMSCEventAction : public G4UserEventAction
{

public:
    BGMSCEventAction();
    ~BGMSCEventAction();

    void BeginOfEventAction(const G4Event* event);
    void EndOfEventAction(const G4Event* event);
};

#endif

    
