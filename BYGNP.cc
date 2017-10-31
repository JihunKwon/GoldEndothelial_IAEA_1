//#undef G4MULTITHREADED

#undef G4VIS_USE

#include <cstdio>
#include <ctime>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4UImanager.hh"
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4ScoringManager.hh"
#include "BGMSCPhysicsList.hh"
#include "BGMSCDetectorConstruction.hh"
#include "BGMSCPrimaryGeneratorAction.hh"
#include "BGMSCEventAction.hh"
#include "BGMSCRunAction.hh"
#include "BGMSCActionInitialization.hh"

#include "G4CsvAnalysisManager.hh"

#include <math.h>

static G4double Sigma;

int main(int argc,char** argv)
{
    // Set the custom seed for the random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4long seed = time(NULL);
    G4Random::setTheSeed(seed);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(1);
#else
    G4RunManager* runManager = new G4RunManager;

//    // Activate UI-command based scorer
//    G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
//    scoringManager->SetVerboseLevel(1);

#endif

    BGMSCDetectorConstruction* massWorld = new BGMSCDetectorConstruction;
    runManager->SetUserInitialization(massWorld);

    G4VModularPhysicsList* physicsList = new BGMSCPhysicsList;
    physicsList->SetVerboseLevel(4);
    runManager->SetUserInitialization(physicsList);

    BGMSCActionInitialization* actionInit = new BGMSCActionInitialization(massWorld);
    runManager->SetUserInitialization(actionInit);
    runManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

//    UImanager->ApplyCommand("/tracking/verbose 1");
//    UImanager->ApplyCommand("/testem/stepMax 5 nm");


#ifdef G4VIS_USE
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  visManager->Initialize();
  UImanager->ApplyCommand("/control/execute init_vis.mac");
  ui->SessionStart();
  delete ui;
  delete visManager;

#endif


  G4UImanager* pUI = G4UImanager::GetUIpointer();

  //pUI->ApplyCommand("/tracking/verbose 1");
  pUI->ApplyCommand("/testem/stepMax 10 nm");

  runManager->BeamOn(50000000); // 25000000 1410065407

  delete runManager;
  return 0;
}
