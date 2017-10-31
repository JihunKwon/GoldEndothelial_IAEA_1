#include "BGMSCPrimaryGeneratorAction.hh"

#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include <G4SingleParticleSource.hh>

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//#define WORLD_SIDE (25. * um)
#define WORLD_SIDE (50. * um)

using namespace CLHEP;

BGMSCPrimaryGeneratorAction::BGMSCPrimaryGeneratorAction()
{
    ////**For Mono Energy**//
//        //Psudo Gamma
//        CircleSourceG = new G4SingleParticleSource;
//        CircleSourceG->SetParticleDefinition(G4Gamma::Definition());//G4Gamma G4Electron


//        // Mono Energy
//        G4SPSEneDistribution* energyG = CircleSourceG->GetEneDist();
//        energyG->SetMonoEnergy(100*keV);
//        // Energy Distribution
////        this->setEnergyToGamma()

//        /////X-Y////Paralell
//        G4SPSPosDistribution* positionG = CircleSourceG->GetPosDist();
//        positionG->SetPosDisType("Plane");
//        positionG->SetPosDisShape("Square");
//        positionG->SetHalfX(WORLD_SIDE/2);
//        positionG->SetHalfY(WORLD_SIDE/2);
//        positionG->SetHalfZ(WORLD_SIDE/2);
//        positionG->SetCentreCoords(G4ThreeVector(0*m, 0*um,  100/2*um));
//        positionG->SetPosRot1(G4ThreeVector(1, 0, 0));
//        positionG->SetPosRot2(G4ThreeVector(0, 1, 0));
//        G4SPSAngDistribution* angularG = CircleSourceG->GetAngDist();
//        angularG->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
//        angularG->DefineAngRefAxes("angref1", G4ThreeVector(0.0,0.0,1.0));
//        CircleSourceG->SetNumberOfParticles(1);

//        /////X-Z////Vertical HPI
//        G4SPSPosDistribution* positionG = CircleSourceG->GetPosDist();
//        positionG->SetPosDisType("Plane");
//        positionG->SetPosDisShape("Square");
//        positionG->SetHalfX(WORLD_SIDE/2);
//        positionG->SetHalfY(WORLD_SIDE/2);
//        positionG->SetHalfZ(25/2*um);
//        positionG->SetCentreCoords(G4ThreeVector(0*m, 100/2*um,  0*um)); // if World is 100um
//        positionG->SetPosRot1(G4ThreeVector(1, 0, 0));
//        positionG->SetPosRot2(G4ThreeVector(0, 0, 1));
//        G4SPSAngDistribution* angularG = CircleSourceG->GetAngDist();
//        angularG->SetParticleMomentumDirection(G4ThreeVector(0.0,-1.0,0.0));
//        angularG->DefineAngRefAxes("angref1", G4ThreeVector(0.0,1.0,0.0));
//        CircleSourceG->SetNumberOfParticles(1);


    ///** For IAEA **///
    G4String fileName = "PSF_-0.05";
    IAEAReader = new G4IAEAphspReader(fileName);

    // Define the rotations
    IAEAReader->SetRotationY(90.*deg); // For parallel irradiation
    //IAEAReader->SetRotationY(0.*deg); // For vertical irradiation

    IAEAReader->SetGlobalPhspTranslation(G4ThreeVector(0,0,100/2*nm)); // 100/2
    IAEAReader->SetTimesRecycled(0);
}

BGMSCPrimaryGeneratorAction::~BGMSCPrimaryGeneratorAction()
{
//    delete CircleSourceG;
    //delete fParticleGun;
    if (IAEAReader)
        delete IAEAReader;
}

void BGMSCPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
//       CircleSourceG->GeneratePrimaryVertex(event);
        //fParticleGun->GeneratePrimaryVertex(event);
       IAEAReader->GeneratePrimaryVertex(event);
}


//G4SPSEneDistribution* BGMSCPrimaryGeneratorAction::setEnergyToBeta() {
//    G4SPSEneDistribution* eneDist = CircleSource->GetEneDist();
//    eneDist->SetEnergyDisType("Arb");
//    eneDist->ArbEnergyHistoFile("Electron_EnergyDistribution_20nm_100keVX-ray_allon_2115_10^6.txt");
//    eneDist->ArbInterpolate("Lin");
//    return eneDist;
//}

//G4SPSEneDistribution* BGMSCPrimaryGeneratorAction::setEnergyToGamma() {
//    G4SPSEneDistribution* eneDistG = CircleSourceG->GetEneDist();
//    eneDistG->SetEnergyDisType("Arb");
//    //eneDistG->ArbEnergyHistoFile("Gamma_EnergyDistribution_allon_10nm_100keVX-ray_2.115_10^9.txt");
//    eneDistG->ArbEnergyHistoFile("Electron_EnergyDistribution_10nm_100keVX-ray.txt");
//    eneDistG->ArbInterpolate("Lin");
//    return eneDistG;
//}




