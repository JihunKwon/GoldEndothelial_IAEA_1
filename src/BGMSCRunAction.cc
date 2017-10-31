#include "BGMSCRunAction.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4TrackingManager.hh"

#include "G4CsvAnalysisManager.hh"

#include <fstream>
#include <sstream>
#include <iostream>

#define _USE_MATH_DEFINES
#define radius 5 // Half of diameter!!

BGMSCRunAction::BGMSCRunAction() : G4UserRunAction()
{}

BGMSCRunAction::~BGMSCRunAction()
{}

void BGMSCRunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void BGMSCRunAction::EndOfRunAction(const G4Run* aRun)
{
}
