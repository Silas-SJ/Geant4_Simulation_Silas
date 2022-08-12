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
/// \file mpdEventAction.cc
/// \brief Implementation of the mpdEventAction class
#include "mpdPrimaryGeneratorAction.hh"
#include "mpdEventAction.hh"
#include "mpdCalorimeterSD.hh"
#include "mpdCalorHit.hh"
#include "mpdDigi.hh"
#include "mpdDigitizer.hh"
#include "mpdAnalysis.hh"
#include "mpdTreeManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4DigiManager.hh"
#include "G4EventManager.hh"
#include "G4VHitsCollection.hh"
#include "Randomize.hh"
#include <iomanip>
#include <set>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdEventAction::mpdEventAction(mpdPrimaryGeneratorAction* ga,mpdTreeManager * tree)
 : G4UserEventAction(),genaction(ga),
   fAbsHCID(-1),
   fGapHCID(-1),
   fScintHCID(-1),
   fDetabsHCID(-1),
   fPmtHCID(-1),
   fPmtHSCID(-1),
   fTreeManager(tree)
  // fcalorSD(0)
{
    mpdDigitizer * mpdDM = new mpdDigitizer( "mpdDigitizer" );
    G4DigiManager::GetDMpointer()->AddNewModule(mpdDM);
   // fcalorSD = new mpdCalorimeterSD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdEventAction::~mpdEventAction()
{

  //delete fcalorSD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorHitsCollection*
mpdEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<mpdCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("mpdEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::EndOfEventAction(const G4Event* event)
{

    
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
    fScintHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("ScintHitsCollection");
    fDetabsHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("DetabsHitsCollection");
    fPmtHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("pmtHitCollection");
    fPmtHSCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("pmtSignalHitCollection");
  }
  
    G4DigiManager * fDM = G4DigiManager::GetDMpointer();
    
  // Get hits collections
  auto absHC = GetHitsCollection(fAbsHCID, event);
  auto gapHC = GetHitsCollection(fGapHCID, event);
  auto scintHC = GetHitsCollection(fScintHCID, event);
  auto detabsHC = GetHitsCollection(fDetabsHCID, event);

  //
  std::map<G4int,G4int> ScintPionDecay;
  std::map<G4int,G4int> ScintMuonDecay;
  std::map<G4int,std::set<G4int>> ScintPionPassed;
  std::map<G4int,std::set<G4int>> ScintMuonPassed; 
  std::map<G4double,G4double> ScintkineticPionEnergy;
  std::map<G4double,G4double> ScintMX_Pion;
  std::map<G4double,G4double> ScintMY_Pion;
  std::map<G4double,G4double> ScintMZ_Pion;
  std::map<G4double,G4double> ScintTheta_Pion;
  std::map<G4double,G4double> ScintPhi_Pion;
  std::set<G4int> ScintPionID;
  std::set<G4int> ScintMuonID;
  std::map<G4int,G4int> Scintn_scint;
  std::map<G4int,G4int> Scintn_cerenk;
  std::map<G4double,G4double> Scintenergy_cerenk;
  std::map<G4double,G4double> Scinttime_cerenk;
  std::map<G4double,G4double> Scintenergy_scint;
  std::map<G4double,G4double> Scinttime_scint;
  //////
  std::map<G4int,G4int> DetabsPionDecay;
  std::map<G4int,G4int> DetabsMuonDecay;
  std::map<G4int,std::set<G4int>> DetabsPionPassed;
  std::map<G4int,std::set<G4int>> DetabsMuonPassed;
  std::set<G4int> DetabsPionID;
  std::set<G4int> DetabsMuonID;
  auto eventID = event->GetEventID();

  for (uint i=0;i<(scintHC->entries());i++){
   auto scintHitlayer = (*scintHC)[i];     

   if(scintHitlayer->GetPionTrackID()!=-111) {
    ScintPionID.insert(scintHitlayer->GetPionTrackID());
    ScintPionPassed[scintHitlayer->GetLayerID()]=ScintPionID;             
   }
   if(scintHitlayer->GetMuonTrackID()!=-111) {
    ScintMuonID.insert(scintHitlayer->GetMuonTrackID());
    ScintMuonPassed[scintHitlayer->GetLayerID()]=ScintMuonID;        
   }
    if(scintHitlayer->GetkineticPionEnergy() != -999999999.){
    ScintkineticPionEnergy[scintHitlayer->GetLayerID()]=scintHitlayer->GetkineticPionEnergy();
  //  std::cout << "Energia passagem1 =" << scintHitlayer->GetkineticPionEnergy() << std::endl;    
   }
   if(scintHitlayer->GetMXPion() != -999999999.){
    ScintMX_Pion[scintHitlayer->GetLayerID()]=scintHitlayer->GetMXPion();
  //  std::cout << "Test passagem1 =" << scintHitlayer->GetMXPion() << std::endl;    
   }
   
   if(scintHitlayer->GetMYPion() != -999999999.){
    ScintMY_Pion[scintHitlayer->GetLayerID()]=scintHitlayer->GetMYPion();
   // std::cout << "Test passagem1 =" << scintHitlayer->GetMYPion() << std::endl;    
   }
   
   if(scintHitlayer->GetMZPion() != -999999999.){
    ScintMZ_Pion[scintHitlayer->GetLayerID()]=scintHitlayer->GetMZPion();
//    std::cout << "Test passagem1 =" << scintHitlayer->GetMZPion() << std::endl;    
   }
   
   if(scintHitlayer->GetThetaPion() != -999999999.){
    ScintTheta_Pion[scintHitlayer->GetLayerID()]=scintHitlayer->GetThetaPion();
   // std::cout << "Test passagem1 =" << scintHitlayer->GetThetaPion() << std::endl;    
   }
   if(scintHitlayer->GetPhiPion() != -999999999.){
    ScintPhi_Pion[scintHitlayer->GetLayerID()]=scintHitlayer->GetPhiPion();
 //   std::cout << "Test passagem1 =" << scintHitlayer->GetPhiPion() << std::endl;    
   }
   ScintPionDecay[scintHitlayer->GetLayerID()]+=scintHitlayer->GetPionDecay();
   ScintMuonDecay[scintHitlayer->GetLayerID()]+=scintHitlayer->GetMuonDecay();
   Scintn_scint[scintHitlayer->GetLayerID()]=scintHitlayer->GetNumberScintillation();
   Scintn_cerenk[scintHitlayer->GetLayerID()]=scintHitlayer->GetNumberCerenkov();
   Scintenergy_cerenk[scintHitlayer->GetLayerID()]=scintHitlayer->GetEnergyCerenkov();
   Scinttime_cerenk[scintHitlayer->GetLayerID()]=scintHitlayer->GetTimeCerenkov();
   Scintenergy_scint[scintHitlayer->GetLayerID()]=scintHitlayer->GetEnergyScintillation();
   Scinttime_scint[scintHitlayer->GetLayerID()]=scintHitlayer->GetTimeScintillation();
  
  }

 

  for (uint i=0;i<(detabsHC->entries());i++){
   auto detabsHitlayer = (*detabsHC)[i];     
   if(detabsHitlayer->GetPionTrackID()!=-111) {
    DetabsPionID.insert(detabsHitlayer->GetPionTrackID());
    DetabsPionPassed[detabsHitlayer->GetLayerID()]=ScintPionID;        
   }
   if(detabsHitlayer->GetMuonTrackID()!=-111) {
    DetabsMuonID.insert(detabsHitlayer->GetMuonTrackID());
    DetabsMuonPassed[detabsHitlayer->GetLayerID()]=DetabsMuonID;        
   }
   DetabsPionDecay[detabsHitlayer->GetLayerID()]+=detabsHitlayer->GetPionDecay();
   DetabsMuonDecay[detabsHitlayer->GetLayerID()]+=detabsHitlayer->GetMuonDecay();
  }
  std::map<G4int,G4int> ScintPionPassou; //PionPassed
  for(auto i=ScintPionPassed.begin();i!=ScintPionPassed.end();++i){
  ScintPionPassou[i->first]=i->second.size();
  // G4cout << "Number of pions: " << i->second.size() << " on detector: " << i->first << G4endl;    
  }
  
  std::map<G4int,G4int> ScintMuonPassou; //MuonPassed
  for(auto i=ScintMuonPassed.begin();i!=ScintMuonPassed.end();++i){
   ScintMuonPassou[i->first]=i->second.size();
  }
  
  std::map<G4int,G4int> DetabsPionPassou; //PionPassed
  for(auto i=DetabsPionPassed.begin();i!=DetabsPionPassed.end();++i){
  DetabsPionPassou[i->first]=i->second.size();
 //  G4cout << "Number of pions: " << i->second.size() << " on detabsorber: " << i->first << G4endl;    
  }
  
  std::map<G4int,G4int> GapPionDecay;
  std::map<G4int,G4int> GapMuonDecay;
  std::map<G4int,std::set<G4int>> GapPionPassed;
  std::map<G4int,std::set<G4int>> GapMuonPassed;
  std::set<G4int> GapPionID;
  std::set<G4int> GapMuonID;
  std::map<G4int,G4int> AbsPionDecay;
  std::map<G4int,G4int> AbsMuonDecay;
  std::map<G4int,std::set<G4int>> AbsPionPassed;
  std::map<G4int,std::set<G4int>> AbsMuonPassed;
  std::set<G4int> AbsPionID;
  std::set<G4int> AbsMuonID;
  std::map<G4double,G4double> AbsXPionPosition;
  std::map<G4double,G4double> AbsYPionPosition;
  std::map<G4double,G4double> AbsZPionPosition;
  for (uint i=0;i<(gapHC->entries());i++){
   auto gapHitlayer = (*gapHC)[i];

   if(gapHitlayer->GetPionTrackID()!=-111) {
    GapPionID.insert(gapHitlayer->GetPionTrackID());
    GapPionPassed[gapHitlayer->GetLayerID()]=GapPionID;
   }
   if(gapHitlayer->GetMuonTrackID()!=-111) {
    GapMuonID.insert(gapHitlayer->GetMuonTrackID());
    GapMuonPassed[gapHitlayer->GetLayerID()]=GapMuonID;
   }
   GapPionDecay[gapHitlayer->GetLayerID()]+=gapHitlayer->GetPionDecay();
   GapMuonDecay[gapHitlayer->GetLayerID()]+=gapHitlayer->GetMuonDecay();
  }

  for (uint i=0;i<(absHC->entries());i++){
   auto absHitlayer = (*absHC)[i];
   if(absHitlayer->GetPionTrackID()!=-111) {
    AbsPionID.insert(absHitlayer->GetPionTrackID());
    AbsPionPassed[absHitlayer->GetLayerID()]=AbsPionID;
   }
   if(absHitlayer->GetMuonTrackID()!=-111) {
    AbsMuonID.insert(absHitlayer->GetMuonTrackID());
    AbsMuonPassed[absHitlayer->GetLayerID()]=AbsMuonID;
   }
   AbsPionDecay[absHitlayer->GetLayerID()]+=absHitlayer->GetPionDecay();
   AbsMuonDecay[absHitlayer->GetLayerID()]+=absHitlayer->GetMuonDecay();
   if(absHitlayer->GetXPionPosition() != -999999999.){
   AbsXPionPosition[absHitlayer->GetLayerID()]=absHitlayer->GetXPionPosition();
   //std::cout << "Test passagem X =" << absHitlayer->GetXPionPosition() << std::endl;    
   }
   if(absHitlayer->GetYPionPosition() != -999999999.){
   AbsYPionPosition[absHitlayer->GetLayerID()]=absHitlayer->GetYPionPosition();
  //  std::cout << "Test passagem Y =" << absHitlayer->GetYPionPosition() << std::endl;    
   }
   if(absHitlayer->GetZPionPosition() != -999999999.){
   AbsZPionPosition[absHitlayer->GetLayerID()]=absHitlayer->GetZPionPosition();
  /// std::cout << "Test passagem Z =" << absHitlayer->GetZPionPosition() << std::endl;    
   }
   
  }
  
  std::map<G4int,G4int> AbsPionPassou; //PionPassed
  for(auto i=AbsPionPassed.begin();i!=AbsPionPassed.end();++i){
   AbsPionPassou[i->first]=i->second.size();
  // G4cout << "Number of pions: " << i->second.size() << " on gap: " << i->first << G4endl;
  }
  
  
  std::map<G4int,G4int> GapPionPassou; //PionPassed
  for(auto i=GapPionPassed.begin();i!=GapPionPassed.end();++i){
   GapPionPassou[i->first]=i->second.size();

  }
  
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

//    PrintEventStatistics(
//      absoHit->GetEdep(), absoHit->GetTrackLength(),
//      gapHit->GetEdep(), gapHit->GetTrackLength());
  }  

 // Digitization
 mpdDigitizer * mpdDM =
 (mpdDigitizer*)fDM->FindDigitizerModule( "mpdDigitizer" );
    if (!mpdDM){
          G4ExceptionDescription msg;
          msg << "Cannot find mpdDigitizer";
          G4Exception("mpdEventAction::EndOfEventAction()",
            "", FatalException, msg);
    }
    
 mpdDM->Digitize();

 G4int DCID = fDM->GetDigiCollectionID("mpdDigitsCollection");
  //  G4cout << "DCID: " << DCID << G4endl;
 mpdDigitsCollection * DC =
 (mpdDigitsCollection*)fDM->GetDigiCollection(DCID);
   // G4cout << "DC->entries() " << DC->entries() << G4endl;

std::vector <G4int> otpv;
std::vector <G4double> ottv;
 if(DC) {
  G4int n_digi =  DC->entries();
  for (G4int i=0;i<n_digi;i++) {

   otpv = (*DC)[i]->GetOverThresholdPMTVec();
   ottv = (*DC)[i]->GetOverThresholdTimeVec();
  }
 }
  // Fill histograms, ntuple
  //

fTreeManager->FillNtuple(
genaction->GetTheta(), genaction->GetPhi(), genaction->GetEnergyPrimary(), genaction->GetMomentumX(), genaction->GetMomentumY(), genaction->GetMomentumZ(), genaction->GetPositionX(), genaction->GetPositionY(), genaction->GetPositionZ(), AbsPionDecay, AbsPionPassou, AbsXPionPosition, AbsYPionPosition, AbsZPionPosition, /*AbsTheta_Pion_Passed, AbsPhi_Pion_Passed,*/ ScintkineticPionEnergy, ScintMX_Pion, ScintMY_Pion, ScintMZ_Pion, ScintTheta_Pion, ScintPhi_Pion, ScintPionDecay, ScintPionPassou, ScintMuonDecay, ScintMuonPassou, Scintn_scint, Scintn_cerenk, Scintenergy_cerenk, Scinttime_cerenk, Scintenergy_scint, Scinttime_scint, DetabsPionDecay, DetabsPionPassou, GapPionDecay, GapPionPassou,  otpv, ottv);


}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
