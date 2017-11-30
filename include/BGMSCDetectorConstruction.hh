#ifndef BGMSCDETECTORCONSTRUCTION_H
#define BGMSCDETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "tls.hh"

class BGMSCDetectorMessenger;
class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////////
struct SCellInfo
{
    int nCellIdx;
};
/////////////////////////////////////////////////////////////////////////////

class BGMSCDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    BGMSCDetectorConstruction();
    ~BGMSCDetectorConstruction() {}

    G4VPhysicalVolume* Construct();

    void SetMaxStep (G4double );
    void SetDistribution(G4String& strDistribution);

 // void ConstructSDandField();
 // This method is used in multi-threaded applications to build
 // per-worker non-shared objects: SensitiveDetectors and Field managers

private:
    G4UserLimits* fStepLimit;   // pointer to user step limits

protected:
    DetectorMessenger* m_pDetectorMessenger;
    G4double m_dWorldSide;
    G4int m_nGnpCount;
    G4String m_strDistribution;
    G4int m_nCellCount;

    void DistributeGnpsSurface (G4LogicalVolume *pWorldLog);
    void DistributeGnpsRandom (G4LogicalVolume *pCubeLog);
    void DistributeGnps4HPI(G4LogicalVolume *pCubeLog);
    void DistributeGnps8HPI(G4LogicalVolume *pCubeLog);
    void DistributeGnps16HPI(G4LogicalVolume *pCubeLog);
    void DistributeGnps24HPI(G4LogicalVolume *pCubeLog);
    void DistributeGnps36HPI(G4LogicalVolume *pCubeLog);
};

#endif

