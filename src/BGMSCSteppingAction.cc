#include "BGMSCSteppingAction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "BGMSCRunAction.hh"

#include "G4CsvAnalysisManager.hh"

#define CellNumb 30
BGMSCSteppingAction::BGMSCSteppingAction()
{
    EdepEndotherialTotal = 0;
    for (int i = 0; i < CellNumb; i++)
    {
        EdepEndotherial[i] = 0;
    }
}

BGMSCSteppingAction::~BGMSCSteppingAction()
{
    //    G4cout << "Total Edep to Total Endotherial cell is " << EdepEndotherialTotal/MeV << G4endl;
        for (int i=0;i<CellNumb;i++)
        {
            G4cout << "EdeptToEndotherial " << i << " " << EdepEndotherial[i] << G4endl;
        }
    //    FILE* fp = fopen("Endotherial_Boxes.csv", "wt");
    //    fprintf(fp, "CellIdx, Edep [MeV]\n");
    //    for (int i=0; i<BoxNumb; i++)
    //    {
    //        fprintf(fp, "%d, %f\n", i, EdepEndotherial[i]);
    //    }
    //    fclose(fp);
    }
    void BGMSCSteppingAction::UserSteppingAction(const G4Step* aStep)
    {

    //    if (aStep->GetTrack()->GetVolume()->GetName() == "EndotherialBoxPhys")
    //    {
    //        EdepEndotherialTotal = EdepEndotherialTotal + aStep->GetTotalEnergyDeposit();
    //    }

    G4TouchableHandle touchable = aStep->GetPreStepPoint()->GetTouchableHandle();		// G4TouchableHandle
    G4VPhysicalVolume* pPhysVolume = touchable->GetVolume();			// G4VPhysicalVolume*

    G4String strPreName = pPhysVolume->GetName();
    if (strPreName == "EndotherialPhys")
    {
        G4int copyNo = pPhysVolume->GetCopyNo();
        EdepEndotherial[copyNo] += aStep->GetTotalEnergyDeposit();
    }
}

