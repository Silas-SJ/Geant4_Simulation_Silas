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
/*
 
 spec-pion.dat
 TF1 * f1 = new TF1("f1","[0]*exp([1]*log([2]*x))",0,4000);
 ****************************************
 Minimizer is Minuit / Migrad
 Chi2                      =  2.80559e-18
 NDf                       =           17
 Edm                       =  2.12624e-26
 NCalls                    =           57
 p0                        =  8.11901e-06   +/-   4.90648e-07
 p1                        =     -1.34047   +/-   0.00919871
 p2                        =     0.378673   +/-   0.016582
 
 */

// 
/// \file mpdPrimaryGeneratorAction.cc
/// \brief Implementation of the mpdPrimaryGeneratorAction class

#include "mpdPrimaryGeneratorAction.hh"
#include "mpdDetectorConstruction.hh"
#include "mpdPrimaryGeneratorMessenger.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath> 
#include <TF1.h>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorAction::mpdPrimaryGeneratorAction(mpdDetectorConstruction* da)
//: G4VUserPrimaryGeneratorAction(), DetConst(da),
 : fParticleGun(nullptr),DetConst(da),// fRndmBeam(0),
 fGunMessenger(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
    //fParticleGun = new G4GeneralParticleSource();
  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun->SetParticleDefinition(particleDefinition);
 // fParticleGun->SetParticleEnergy(0.001*GeV);
  //fParticleGun->SetParticlePosition(G4ThreeVector(100., 30., -38.));
  fGunMessenger = new mpdPrimaryGeneratorMessenger(this);  
// Create muon (pion) angular spectrum //
  for (int n = 0; n < 150; n++)
  { 
    thetabin[n] = cos((M_PI/2.*n)/200.)*cos((M_PI/2.*n)/200.);
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorAction::~mpdPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }
  
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("mpdPrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 


    CLHEP::RandGeneral muonDist(thetabin,150);
    Theta = muonDist.shoot()*M_PI/2;
   // Theta = 0.;
    std::cout<< "Theta ="<< Theta << std::endl;
    G4double costheta = cos(Theta); 
    G4double sintheta = std::sqrt(1. - pow(costheta,2));
    Phi = CLHEP::twopi*G4UniformRand();
   //Phi = 0.;
    std::cout<< "Phi ="<< Phi << std::endl;
    pX = sintheta*cos(Phi); // 
    pY = sintheta*sin(Phi);
    pZ = costheta;
    
   
 //#------------- Above Construction------------#
  X =-110.*m -45*(2*G4UniformRand()-1)*m;
  Y =-10.5*m -45*(2*G4UniformRand()-1)*m;
  Z = -22.2*m; //(above construction)
  fParticleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
  
  std::cout<< "XPosition ="<< X << "mm"<< std::endl;
  std::cout<< "YPosition ="<< Y << "mm"<< std::endl;
  std::cout<< "ZPosition ="<< Z << "mm"<< std::endl;
  
  
 // # :::::: spectrum energy fit from spec-pion.dat ::::://
  TF1 * f1 = new TF1("f1","[0]*pow(x,[1])",0.1,2);
  f1->SetParameter(0,5.68934e-5);
  f1->SetParameter(1,-2.56427);
  double Energy = f1->GetRandom()*GeV;
  fParticleGun->SetParticleEnergy(Energy);
  
  if (X>=-65*m && X<=-115*m && Y>=12.5*m && Y<=-12.5*m){
      return;
  }
  
  else if (X<=-114.5*m && pX<0){
      return ; 
  }
  else if (X>=-64.5*m && pX >0){
     return ;
  }
  else if (Y<=-12*m && pY <0 ){
     return ;
  }
  else if (Y>=12*m && pY >0){
      return ;
  }
  else{
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pX, pY, pZ));
    }
 
std::cout << "Energy Started = " << Energy << " MeV = "<<  std::endl;
 
/*
//#------------- Above Detector------------#
  X =-110.*m -0.20*(2*G4UniformRand()-1)*m;
  Y = -10.5*m -0.20*(2*G4UniformRand()-1)*m;
  Z = 18.5*m/2 + 0.77*m;  //(above first scintillator) +0,78 is inside scintillator //
 // Z = -22.2*m; //(above construction)

 
  fParticleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
 
  if(Y<=-10.699*m && pY<0 ){
   return ; 
  }
  else if(X>=-109.799*m && pX >0){ 
   return ; 
  }
  else if(Y>=-10.299*m && pY>0 ){
   return ; 
  }
   else if(X<=-110.199*m && pX<0){
   return ; 
  }
  else{
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pX, pY, pZ));
  }

  // :::::::Randomizando a energia:::::::::: //
   G4double E0 = 10 *MeV;
   G4double Ef = 100 *MeV;
   G4double energystart = E0 + G4UniformRand() * (Ef - E0) ;
   fParticleGun->SetParticleEnergy(energystart); 
   std::cout << "Gun Energy =" << energystart << std::endl; 
*/



   gunenergy = fParticleGun->GetParticleEnergy();/// keep at end

  

  fParticleGun->GeneratePrimaryVertex(anEvent); // keep at end
}
void mpdPrimaryGeneratorAction::SetTheta (G4double val ) 
   { Theta = val;}
  
   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

