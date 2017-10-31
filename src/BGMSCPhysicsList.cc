#include "BGMSCPhysicsList.hh"

// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4ParallelWorldScoringProcess.hh"
#include "G4EmParameters.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4WentzelVIModel.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4UAtomicDeexcitation.hh"

#include "StepMax.hh"

using namespace CLHEP;

BGMSCPhysicsList::BGMSCPhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1*nm; //0.1*nm
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  stepMaxProcess  = 0;

  G4EmParameters* emParameters = G4EmParameters::Instance();
  emParameters->SetMinEnergy(0*eV);
  emParameters->SetMaxEnergy(100*MeV);
  emParameters->SetNumberOfBinsPerDecade(50);

//  emParameters->SetMuHadLateralDisplacement(true);
//  emParameters->SetLateralDisplacement(true);

//  emParameters->SetMscStepLimitType(fUseDistanceToBoundary);
//  emParameters->SetMscGeomFactor(2.5);
//  emParameters->SetMscRangeFactor(0.04);

  emParameters->SetApplyCuts(true);
  emParameters->SetVerbose(true);

  SetVerboseLevel(1);

  //RegisterPhysics(new G4StepLimit)
  RegisterPhysics(new G4EmPenelopePhysics);
  RegisterPhysics(new G4StepLimiterPhysics);
  // RegisterPhysics(new PhysListEmStandardSingleSc);
  RegisterPhysics(new G4HadronPhysicsQGSP_BIC);
  RegisterPhysics(new G4EmExtraPhysics);
  RegisterPhysics(new G4HadronElasticPhysics);
  RegisterPhysics(new G4StoppingPhysics);
  RegisterPhysics(new G4IonBinaryCascadePhysics);
  RegisterPhysics(new G4NeutronTrackingCut);
}

BGMSCPhysicsList::~BGMSCPhysicsList()
{}

void BGMSCPhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    //Deexcitation
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    de->SetFluo(true);
    de->SetAuger(true);
    de->SetPIXE(true);
    de->SetAugerCascade(true);
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);

    // step limitation (as a full process)
    //
    AddStepMax();
}

//void PhysicsLists::AddP

void BGMSCPhysicsList::SetCuts()
{
 // set cut values for gamma at first and for e- second and next for e+,
 // because some processes for e+/e- need cut values for gamma
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);
}

void BGMSCPhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    stepMaxProcess = new StepMax();

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while ((*theParticleIterator)()){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();

        if (stepMaxProcess->IsApplicable(*particle) && pmanager)
        {
            pmanager->AddDiscreteProcess(stepMaxProcess);
        }
    }
}
