#include "BGMSCDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4VSolid.hh"
#include "G4Sphere.hh"
#include "Randomize.hh"
#include "G4UserLimits.hh"

#include "math.h"
using namespace CLHEP;

/*
20nm
10mg 388199
15mg 582298
20mg 776398
30mg 1164596
40mg 1552795
 *
100nm
10mg 3106
20mg 6211
30mg 9317
40mg 12422
*/
///////////////////////////////////////////////////////////////////
#define WORLD_SIDE (100. * um) // 100 for 8HPI
//#define WORLD_SIDE (25. * um) for no HPI study
#define VESSEL_OUTER_DIAM (20.0 * um)
#define VESSEL_INNER_DIAM (16.0 * um)
#define VESSEL_HIGHT (10.0 * um)
#define GNP_DIAM (100 * nm) // diameter!!
#define GNP_COUNT 12422// Number! 15528*0.7176=11143
#define CELL_COUNT 30
#define SHELL1_DIAM (20.0 * um)
#define SHELL2_DIAM (40.0 * um)
#define SHELL3_DIAM (60.0 * um)
#define SHELL4_DIAM (80.0 * um)
#define SHELL5_DIAM (100.0 * um)


//4HPI
#define SHELL_FRC20NM_4HPI_0_10 0.6345
#define SHELL_FRC20NM_4HPI_10_20 0.2886
#define SHELL_FRC20NM_4HPI_20_30 0.0657
#define SHELL_FRC20NM_4HPI_30_40 0.0100
#define SHELL_FRC20NM_4HPI_40_50 0.0011
#define SHELL_FRC100NM_4HPI_0_10 0.7180
#define SHELL_FRC100NM_4HPI_10_20 0.2379
#define SHELL_FRC100NM_4HPI_20_30 0.0394
#define SHELL_FRC100NM_4HPI_30_40 0.0044
#define SHELL_FRC100NM_4HPI_40_50 0.0004

//8HPI
#define SHELL_FRC20NM_8HPI_0_10 0.4026
#define SHELL_FRC20NM_8HPI_10_20 0.2517
#define SHELL_FRC20NM_8HPI_20_30 0.1507
#define SHELL_FRC20NM_8HPI_30_40 0.1037
#define SHELL_FRC20NM_8HPI_40_50 0.0913
#define SHELL_FRC100NM_8HPI_0_10 0.5155
#define SHELL_FRC100NM_8HPI_10_20 0.3039
#define SHELL_FRC100NM_8HPI_20_30 0.0803
#define SHELL_FRC100NM_8HPI_30_40 0.0539
#define SHELL_FRC100NM_8HPI_40_50 0.0464

//16HPI
#define SHELL_FRC20NM_16HPI_0_10 0.1621
#define SHELL_FRC20NM_16HPI_10_20 0.2949
#define SHELL_FRC20NM_16HPI_20_30 0.2683
#define SHELL_FRC20NM_16HPI_30_40 0.1628
#define SHELL_FRC20NM_16HPI_40_50 0.0740
#define SHELL_FRC100NM_16HPI_0_10 0.2658
#define SHELL_FRC100NM_16HPI_10_20 0.3522
#define SHELL_FRC100NM_16HPI_20_30 0.2333
#define SHELL_FRC100NM_16HPI_30_40 0.1031
#define SHELL_FRC100NM_16HPI_40_50 0.0341

//24HPI
#define SHELL_FRC20NM_24HPI_0_10 0.0652
#define SHELL_FRC20NM_24HPI_10_20 0.1781
#define SHELL_FRC20NM_24HPI_20_30 0.2431
#define SHELL_FRC20NM_24HPI_30_40 0.2212
#define SHELL_FRC20NM_24HPI_40_50 0.1509
#define SHELL_FRC100NM_24HPI_0_10 0.1370
#define SHELL_FRC100NM_24HPI_10_20 0.2723
#define SHELL_FRC100NM_24HPI_20_30 0.2707
#define SHELL_FRC100NM_24HPI_30_40 0.1793
#define SHELL_FRC100NM_24HPI_40_50 0.0891

//36HPI
#define SHELL_FRC20NM_36HPI_0_10 0.0167
#define SHELL_FRC20NM_36HPI_10_20 0.0682
#define SHELL_FRC20NM_36HPI_20_30 0.1397
#define SHELL_FRC20NM_36HPI_30_40 0.1907
#define SHELL_FRC20NM_36HPI_40_50 0.1952
#define SHELL_FRC100NM_36HPI_0_10 0.0507
#define SHELL_FRC100NM_36HPI_10_20 0.1512
#define SHELL_FRC100NM_36HPI_20_30 0.2254
#define SHELL_FRC100NM_36HPI_30_40 0.2240
#define SHELL_FRC100NM_36HPI_40_50 0.1670

///////////////////////////////////////////////////////////////////

BGMSCDetectorConstruction::BGMSCDetectorConstruction()
    :fStepLimit(NULL)
{
    m_dWorldSide = WORLD_SIDE;
    m_nGnpCount = GNP_COUNT;
    m_strDistribution = "36HPI";   // Constrained  Random None 8HPI 24HPI
    m_nCellCount = CELL_COUNT;
}

G4VPhysicalVolume* BGMSCDetectorConstruction::Construct()
{
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    G4VisAttributes* pVisAttributesEndotherial = new G4VisAttributes;
    pVisAttributesEndotherial->SetForceWireframe(true);
    pVisAttributesEndotherial->SetForceAuxEdgeVisible(true);
    pVisAttributesEndotherial->SetForceSolid(true);
    pVisAttributesEndotherial->SetVisibility(true);
    pVisAttributesEndotherial->SetColor(1.0, 0.0, 0.0); // red

//    G4VisAttributes* pVisAttributesShell = new G4VisAttributes;
//    pVisAttributesShell->SetForceWireframe(true);
//    pVisAttributesShell->SetForceAuxEdgeVisible(true);
//    pVisAttributesShell->SetForceSolid(false);
//    pVisAttributesShell->SetVisibility(true);
//    pVisAttributesShell->SetColor(1.0, 1.0, 1.0); // cyan


    // Build materials
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* Au = nistManager->FindOrBuildMaterial("G4_Au");
    G4Material* vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
    G4Material* water = nistManager->FindOrBuildMaterial("G4_WATER");

    // Define World
    G4Box* pWorldBox = new G4Box("WorldBox", WORLD_SIDE/2, WORLD_SIDE/2, WORLD_SIDE/2);
    G4LogicalVolume *pWorldLogic = new G4LogicalVolume(pWorldBox, water, "WorldLog");
    G4VPhysicalVolume *pWorldPhys = new G4PVPlacement(0, G4ThreeVector(), pWorldLogic, "WorldPhys", 0, false, 0);

//    // Define Endotherial cell
//    G4Tubs* pEndotherialTubs = new G4Tubs("EndotherialTubs", VESSEL_INNER_DIAM/2, VESSEL_OUTER_DIAM/2, VESSEL_HIGHT/2, 0*deg, 360*deg);
//    G4LogicalVolume *pEndotherialLogic = new G4LogicalVolume(pEndotherialTubs, water, "EndotherialLog");
//    G4VPhysicalVolume *pEndotherialPhys = new G4PVPlacement(0, G4ThreeVector(), pEndotherialLogic, "EndotherialPhys", pWorldLogic, 0, 0);
//    pEndotherialLogic->SetVisAttributes(pVisAttributesEndotherial);

//    // Psudo Shells
//    for (int ShellNo = 1; ShellNo <= 5; ShellNo++)
//    {
//        G4Tubs* pShellTubs = new G4Tubs("ShellTubes", 0, VESSEL_OUTER_DIAM/2*ShellNo, VESSEL_HIGHT/2, 0*deg, 360*deg);
//        G4LogicalVolume *pShellLogic = new G4LogicalVolume(pShellTubs, water, "ShellLog");
//        G4VPhysicalVolume *pShellPhys = new G4PVPlacement(0, G4ThreeVector(), pShellLogic, "ShellPhys", pWorldLogic, 0, 0);
//        pShellLogic->SetVisAttributes(pVisAttributesShell);
//    }

    // Parameterize each Endotherial cell.
    for (int nCellIdx = 0; nCellIdx < m_nCellCount; nCellIdx++)
    {
        G4RotationMatrix* rotm = new G4RotationMatrix;
        rotm->rotateZ((nCellIdx+1) * 12*deg - 90*deg); //Set the first position 0 oclock., clock wise.
                                                                                                                               /* Cell bin size */
        G4Tubs* pEndotherialTubs = new G4Tubs("EndotherialTubs", VESSEL_INNER_DIAM/2, VESSEL_OUTER_DIAM/2, VESSEL_HIGHT/2, 0*deg, 12*deg); // Last deg is the size of cell!
        G4LogicalVolume *pEndotherialLogic = new G4LogicalVolume(pEndotherialTubs, water, "EndotherialLog");
        G4VPhysicalVolume *EndotherialPhys = new G4PVPlacement(rotm, G4ThreeVector(), pEndotherialLogic, "EndotherialPhys", pWorldLogic, 0, nCellIdx);
        pEndotherialLogic->SetVisAttributes(pVisAttributesEndotherial);
    }

    if (m_strDistribution == "Constrained")
    {
        DistributeGnpsSurface(pWorldLogic);
    }
    else if (m_strDistribution == "Random")
    {
        DistributeGnpsRandom(pWorldLogic);
    }
    else if (m_strDistribution == "4HPI")
    {
        DistributeGnps4HPI(pWorldLogic);
    }
    else if (m_strDistribution == "8HPI")
    {
        DistributeGnps8HPI(pWorldLogic);
    }
    else if (m_strDistribution == "16HPI")
    {
        DistributeGnps16HPI(pWorldLogic);
    }
    else if (m_strDistribution == "24HPI")
    {
        DistributeGnps24HPI(pWorldLogic);
    }
    else if (m_strDistribution == "36HPI")
    {
        DistributeGnps36HPI(pWorldLogic);
    }
    else if (m_strDistribution == "None")
    {
    }
    else
    {
        G4Exception("Invalid distribution name", "001", G4ExceptionSeverity::FatalException, "Comment");
    }

    G4double maxStep = 1*nm;
    fStepLimit = new G4UserLimits(maxStep);

    return pWorldPhys;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
///////////////////////////CONSTRAINED////////////////////////////////////////
// Distribute GNP randomly on the inner surface of the blood vessel.
void BGMSCDetectorConstruction::DistributeGnpsSurface(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributes = new G4VisAttributes;
    pVisAttributes->SetForceWireframe(true);
    pVisAttributes->SetForceAuxEdgeVisible(true);
    pVisAttributes->SetForceSolid(false);
    pVisAttributes->SetVisibility(true);
    pVisAttributes->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
    {
        retry:
        // Compute a random position for the GNP
        double dTheta = G4UniformRand() * 360.*deg;
        double dGnpX = ((VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
        double dGnpY = ((VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
        double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Special case for one GNP at the center.
                if (m_nGnpCount == 1)
                {
                    dGnpX = 0.0;
                    dGnpY = 0.0;
                    dGnpZ = 0.0;
                }

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry;
                }
                printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

                G4double r = VESSEL_INNER_DIAM/2.;
                G4double Cos = cos(dGnpX/r);
                G4double theta = acos(Cos)*180/M_PI;
//                if (dGnpX<0 && dGnpY>0) theta = theta+90; //2nd
//                else if(dGnpX<0 && dGnpY<0) theta = theta+180; //3rd
//                else if(dGnpX>0 && dGnpY<0) theta = theta+270; //4th

                for (int i=0; i<60; i++)
                {
                    if((i<=theta) && (theta<(i+1)))
                    {
                        ary[i] += 1;
                    }
                }

                int nCopyNumber = nGnpIdx; //Why?

                new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);

                pGnpLog->SetVisAttributes(pVisAttributes);

                aryGnpInfo[nGnpIdx].dPosX = dGnpX;
                aryGnpInfo[nGnpIdx].dPosY = dGnpY;
                aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
    }
    FILE* fp =fopen("position.txt", "wt");
    for (int i=0; i<60; i++)
    {
        printf("%d %d\n", i, ary[i]);
        fprintf(fp, "%d %d\n", i, ary[i]);
    }
    fclose(fp);
    delete [] aryGnpInfo;
}

///////////////////////////RANDOM////////////////////////////////////////
// Distribute GNP randomly in the blood vessel.
void BGMSCDetectorConstruction::DistributeGnpsRandom(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributes = new G4VisAttributes;
    pVisAttributes->SetForceWireframe(true);
    pVisAttributes->SetForceAuxEdgeVisible(true);
    pVisAttributes->SetForceSolid(false);
    pVisAttributes->SetVisibility(true);
    pVisAttributes->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...");

    for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
    {
        retry:
        // Compute a random position for the GNP
        double dTheta = G4UniformRand() * 2*M_PI;
        double dRand = G4UniformRand();
        double dGnpX = (sqrt(dRand) * (VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
        double dGnpY = (sqrt(dRand) * (VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
        double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Special case for one GNP at the center.
                if (m_nGnpCount == 1)
                {
                    dGnpX = 0.0;
                    dGnpY = 0.0;
                    dGnpZ = 0.0;
                }

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
//                    if (fabs(dx) > GNP_DIAM) continue;
//                    if (fabs(dy) > GNP_DIAM) continue;
//                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry;
                }
                printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

                int nCopyNumber = nGnpIdx; //Why?

                new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);

                pGnpLog->SetVisAttributes(pVisAttributes);

                aryGnpInfo[nGnpIdx].dPosX = dGnpX;
                aryGnpInfo[nGnpIdx].dPosY = dGnpY;
                aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
    }
    delete [] aryGnpInfo;
}

///////////////////////////4HPI////////////////////////////////////////
// Distribute GNP with migration after 4 hour post injection (8HPI)
void BGMSCDetectorConstruction::DistributeGnps4HPI(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
    pVisAttributesGold->SetForceWireframe(true);
    pVisAttributesGold->SetForceAuxEdgeVisible(true);
    pVisAttributesGold->SetForceSolid(false);
    pVisAttributesGold->SetVisibility(true);
    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    if (GNP_DIAM == 20 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_4HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC20NM_4HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC20NM_4HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC20NM_4HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC20NM_4HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry20_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry20_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry20_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry20_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry20_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry20_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry20_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry20_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry20_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
    if (GNP_DIAM == 100 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC100NM_4HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC100NM_4HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC100NM_4HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC100NM_4HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC100NM_4HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry100_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry100_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry100_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry100_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry100_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry100_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry100_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry100_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry100_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
//    FILE* fp =fopen("position.txt", "wt");
//    for (int i=0; i<60; i++)
//    {
//        printf("%d %d\n", i, ary[i]);
//        fprintf(fp, "%d %d\n", i, ary[i]);
//    }
//    fclose(fp);
    delete [] aryGnpInfo;
}

///////////////////////////8HPI////////////////////////////////////////
// Distribute GNP with migration after 8 hour post injection (8HPI)
void BGMSCDetectorConstruction::DistributeGnps8HPI(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
    pVisAttributesGold->SetForceWireframe(true);
    pVisAttributesGold->SetForceAuxEdgeVisible(true);
    pVisAttributesGold->SetForceSolid(false);
    pVisAttributesGold->SetVisibility(true);
    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    if (GNP_DIAM == 20 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_8HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC20NM_8HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC20NM_8HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC20NM_8HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC20NM_8HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry20_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry20_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry20_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry20_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry20_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry20_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry20_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry20_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry20_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
    if (GNP_DIAM == 100 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC100NM_8HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC100NM_8HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC100NM_8HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC100NM_8HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC100NM_8HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry100_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry100_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry100_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry100_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry100_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry100_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry100_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry100_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry100_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
//    FILE* fp =fopen("position.txt", "wt");
//    for (int i=0; i<60; i++)
//    {
//        printf("%d %d\n", i, ary[i]);
//        fprintf(fp, "%d %d\n", i, ary[i]);
//    }
//    fclose(fp);
    delete [] aryGnpInfo;
}


///////////////////////////16HPI////////////////////////////////////////
// Distribute GNP with migration after 8 hour post injection (16HPI)
void BGMSCDetectorConstruction::DistributeGnps16HPI(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
    pVisAttributesGold->SetForceWireframe(true);
    pVisAttributesGold->SetForceAuxEdgeVisible(true);
    pVisAttributesGold->SetForceSolid(false);
    pVisAttributesGold->SetVisibility(true);
    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    if (GNP_DIAM == 20 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_16HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC20NM_16HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC20NM_16HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC20NM_16HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC20NM_16HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry20_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry20_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry20_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry20_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry20_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry20_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry20_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry20_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry20_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
    if (GNP_DIAM == 100 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC100NM_16HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC100NM_16HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC100NM_16HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC100NM_16HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC100NM_16HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry100_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry100_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry100_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry100_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry100_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry100_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry100_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry100_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry100_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
//    FILE* fp =fopen("position.txt", "wt");
//    for (int i=0; i<60; i++)
//    {
//        printf("%d %d\n", i, ary[i]);
//        fprintf(fp, "%d %d\n", i, ary[i]);
//    }
//    fclose(fp);
    delete [] aryGnpInfo;
}


///////////////////////////24HPI////////////////////////////////////////
// Distribute GNP with migration after 24 hour post injection (24HPI)
void BGMSCDetectorConstruction::DistributeGnps24HPI(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

   // m_nGnpCount = GNP_COUNT*(SHELL_FRC100NM_24HPI_0_10+SHELL_FRC100NM_24HPI_10_20+SHELL_FRC100NM_24HPI_20_30
   //                          +SHELL_FRC100NM_24HPI_30_40+SHELL_FRC100NM_24HPI_40_50);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
    pVisAttributesGold->SetForceWireframe(true);
    pVisAttributesGold->SetForceAuxEdgeVisible(true);
    pVisAttributesGold->SetForceSolid(false);
    pVisAttributesGold->SetVisibility(true);
    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    if (GNP_DIAM == 20 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_24HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC20NM_24HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC20NM_24HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC20NM_24HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC20NM_24HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry20_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry20_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry20_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry20_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry20_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry20_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry20_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry20_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry20_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }

    if (GNP_DIAM == 100 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC100NM_24HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC100NM_24HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC100NM_24HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC100NM_24HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC100NM_24HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry100_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry100_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry100_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry100_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry100_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry100_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry100_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry100_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry100_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
//    FILE* fp =fopen("position.txt", "wt");
//    for (int i=0; i<60; i++)
//    {
//        printf("%d %d\n", i, ary[i]);
//        fprintf(fp, "%d %d\n", i, ary[i]);
//    }
//    fclose(fp);
    delete [] aryGnpInfo;
}


///////////////////////////36HPI////////////////////////////////////////
// Distribute GNP with migration after 36 hour post injection (36HPI)
void BGMSCDetectorConstruction::DistributeGnps36HPI(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
    pVisAttributesGold->SetForceWireframe(true);
    pVisAttributesGold->SetForceAuxEdgeVisible(true);
    pVisAttributesGold->SetForceSolid(false);
    pVisAttributesGold->SetVisibility(true);
    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...\n");
    G4int ary[360]={0};

    if (GNP_DIAM == 20 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_36HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC20NM_36HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC20NM_36HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC20NM_36HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC20NM_36HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry20_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry20_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry20_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry20_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry20_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry20_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry20_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry20_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry20_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry20_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }

    if (GNP_DIAM == 100 * nm)
    {
        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
        {
            double dTheta, dRand, dGnpX, dGnpY, dGnpZ; //Declear to use at the end "fprintf()"
            double dGnpR; // Radial Direction

            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC100NM_36HPI_0_10;
            G4int GNP_COUNT_10_20 = GNP_COUNT*SHELL_FRC100NM_36HPI_10_20;
            G4int GNP_COUNT_20_30 = GNP_COUNT*SHELL_FRC100NM_36HPI_20_30;
            G4int GNP_COUNT_30_40 = GNP_COUNT*SHELL_FRC100NM_36HPI_30_40;
            G4int GNP_COUNT_40_50 = GNP_COUNT*SHELL_FRC100NM_36HPI_40_50;

            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
            {
                retry100_1:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL1_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_1;
                }
            }
            else if (GNP_COUNT_0_10 <= nGnpIdx && nGnpIdx < (GNP_COUNT_0_10 + GNP_COUNT_10_20)) // Shell2
            {
                retry100_2:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL2_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                dGnpR = sqrt(dGnpX*dGnpX + dGnpY*dGnpY)/um;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL1_DIAM/2.)
                    goto retry100_2;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_2;
                }
            }
            else if ((GNP_COUNT_0_10 + GNP_COUNT_10_20) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30)) // Shell3
            {
                retry100_3:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL3_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL3_DIAM)
                    goto retry100_3;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_3;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40)) // Shell4
            {
                retry100_4:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL4_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL4_DIAM)
                    goto retry100_4;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_4;
                }
            }
            else if ((GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40) <= nGnpIdx &&
                     nGnpIdx < (GNP_COUNT_0_10+GNP_COUNT_10_20+GNP_COUNT_20_30+GNP_COUNT_30_40+GNP_COUNT_40_50)) // Shell5
            {
                retry100_5:
                // Compute a random position for the GNP
                dTheta = G4UniformRand() * 2*M_PI;
                dRand = G4UniformRand();
                dGnpX = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
                dGnpY = (sqrt(dRand) * (SHELL5_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
                dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;
                // Check if this GNP is in the inner shell
                if (dGnpR <= SHELL5_DIAM)
                    goto retry100_5;

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry100_5;
                }
            }
            else { G4cout << "Shell index error!!" << G4endl;}

            printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

            int nCopyNumber = nGnpIdx;
            new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
            pGnpLog->SetVisAttributes(pVisAttributesGold);

            aryGnpInfo[nGnpIdx].dPosX = dGnpX;
            aryGnpInfo[nGnpIdx].dPosY = dGnpY;
            aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
        }
    }
//    FILE* fp =fopen("position.txt", "wt");
//    for (int i=0; i<60; i++)
//    {
//        printf("%d %d\n", i, ary[i]);
//        fprintf(fp, "%d %d\n", i, ary[i]);
//    }
//    fclose(fp);
    delete [] aryGnpInfo;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BGMSCDetectorConstruction::SetMaxStep(G4double maxStep)
{
    if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}
