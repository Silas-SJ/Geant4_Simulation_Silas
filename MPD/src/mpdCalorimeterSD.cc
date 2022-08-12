//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file mpdCalorimeterSD.cc
/// \brief Implementation of the mpdCalorimeterSD class

#include "mpdCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include <vector>
#include <cmath>
#include "G4SystemOfUnits.hh" 
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorimeterSD::mpdCalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells
                            )
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorimeterSD::~mpdCalorimeterSD()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new mpdCalorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

G4RunManager *rm = G4RunManager::GetRunManager();
eventID= rm->GetCurrentEvent()->GetEventID();
oldID = -9999;
kineticPionEnergy = -999999999.;
MX_pion = -999999999.;
MY_pion = -999999999.;
MZ_pion = -999999999.;
Theta_pion = -999999999.;
Phi_pion = -999999999.;
M_Pion_total = -999999999.;
X_pion= -999999999.;
Y_pion= -999999999.;
Z_pion= -999999999.;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool mpdCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
  // energy deposit
   
  auto edep = step->GetTotalEnergyDeposit();
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0. ) return false;
    stepLength = step->GetStepLength();

 // if ( edep==0. && stepLength == 0. ) return false;

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id 
  auto layerNumber = touchable->GetReplicaNumber(1);
  
   mpdCalorHit * hit = new mpdCalorHit(); 
    if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("mpdCalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
    }

  // Get hit for total accounting
//  auto hitTotal 
//    = (*fHitsCollection)[fHitsCollection->entries()-1];
  auto piondecay = false;
  auto muondecay = false;
  auto pionpassed = false;
  auto muonpassed = false;
  auto insidVolume = false;
  //G4int detNo = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
 
  ParentID = step->GetTrack()->GetParentID();
  
  if(abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211) {
  
  if (step->GetTrack()->GetTrackID() != oldID && ParentID != oldID) {
     pionpassed = true;

    if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Scinti"){
      insidVolume = true;

      kineticPionEnergy = step->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/MeV;  
      
      MX_pion = step->GetTrack()->GetMomentumDirection().x();

      MY_pion = step->GetTrack()->GetMomentumDirection().y();

      
      MZ_pion = step->GetTrack()->GetMomentumDirection().z();

      
      M_Pion_total = std::sqrt(pow(MX_pion,2) + pow(MY_pion,2) + pow(MZ_pion,2));
      
      Theta_pion = (std::acos(MZ_pion/M_Pion_total))*(180/(CLHEP::twopi/2)); ///It is in deg

      
      Phi_pion = (std::atan(MY_pion/MX_pion))*(180/(CLHEP::twopi/2)); ///It is in deg

   }
  } 
   oldID = step->GetTrack()->GetTrackID();
  
   }
   
   

   
  if(abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 13){
   muonpassed = true;

  }
  
 if ((step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Scintillation"||step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay")&& abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211){  
      piondecay=true;

      X_pion= (step->GetPostStepPoint()->GetPosition().getX())/m;

      Y_pion= (step->GetPostStepPoint()->GetPosition().getY())/m;

      Z_pion= (step->GetPostStepPoint()->GetPosition().getZ())/m;
  
    }
   
  if ((step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Scintillation"||step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay") && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 13){
      muondecay=true;
      hit->AddMuonDecay();
    }


   auto proc_man =
      step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetProcessManager();
    G4ProcessVector* proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);
    G4int n_proc              = proc_vec->entries();

    G4int n_scint = 0;
    G4int n_cerenk   = 0;
    
   if(abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211){
    for(G4int i = 0; i < n_proc; ++i)
    {
      G4String proc_name = (*proc_vec)[i]->GetProcessName();
      if(proc_name.compare("Cerenkov") == 0)
      {
        auto cerenk = (G4Cerenkov*) (*proc_vec)[i];
        n_cerenk = cerenk->GetNumPhotons();
      }
      else if(proc_name.compare("Scintillation") == 0)
      {
        auto scint = (G4Scintillation*) (*proc_vec)[i];
        n_scint  = (scint->GetNumPhotons());
      }
   }

    }
     // loop over secondaries, create statistics
    const std::vector<const G4Track*>* secondaries =
      step->GetSecondaryInCurrentStep();
      
    G4double energy_cerenk = 0.;
    G4double time_cerenk = 0.;
    G4double energy_scint = 0.;
    G4double time_scint = 0.;
    
    for(auto sec : *secondaries)
    {
      if(sec->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) // G4OpticalPhoton) //opticalphoton
      {
        G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
        if(creator_process.compare("Cerenkov") == 0)
        {
          energy_cerenk = (sec->GetKineticEnergy()) / eV;
          time_cerenk = (sec->GetGlobalTime()) / ns;
        }
        else if(creator_process.compare("Scintillation") == 0)
        {
          energy_scint = (sec->GetKineticEnergy()) / eV;
          time_scint = (sec->GetGlobalTime()) / ns;

        }
      }
    }
 
   hit->Add(edep, stepLength, layerNumber, n_scint, n_cerenk, energy_cerenk, time_cerenk, energy_scint, time_scint);
   
    if(pionpassed){
      hit->SetPionTrackID(step->GetTrack()->GetTrackID());
    }
    if(muonpassed){
      hit->SetMuonTrackID(step->GetTrack()->GetTrackID());
    }
    if(piondecay){
      hit->AddPionDecay();
      hit->SetX_PionPosition(X_pion); 
      hit->SetY_PionPosition(Y_pion);
      hit->SetZ_PionPosition(Z_pion);
    }
    if(muondecay){
    hit->AddMuonDecay();
    }
    if(insidVolume){
      hit->SetkineticPionEnergy(kineticPionEnergy);
      hit->SetMXPion(MX_pion); 
      hit->SetMYPion(MY_pion);
      hit->SetMZPion(MZ_pion);
      hit->SetThetaPion(Theta_pion);
      hit->SetPhiPion(Phi_pion);
    }
    fHitsCollection->insert(hit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     auto nofHits = fHitsCollection->entries();
    /* 
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
       */
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
