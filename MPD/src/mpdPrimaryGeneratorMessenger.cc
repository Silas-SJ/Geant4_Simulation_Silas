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
/// \file electromagnetic/TestEm12/src/PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "mpdPrimaryGeneratorMessenger.hh"

#include "mpdPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#define CHECKPG()

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorMessenger::mpdPrimaryGeneratorMessenger(
                                                   mpdPrimaryGeneratorAction* Gun)
:G4UImessenger(),fAction(Gun),
 fGunDir(0),  
 thetaCmd(0)
 

{ 
  fGunDir = new G4UIdirectory("/gun2/");
  fGunDir->SetGuidance("gun control");
  //# Theta ang primary particle #
  thetaCmd = new G4UIcmdWithADoubleAndUnit("/gun2/theta",this);
  thetaCmd->SetGuidance("Set theta");
  thetaCmd->SetParameterName("theta",false);
  thetaCmd->SetDefaultUnit("deg");
  // thetaCmd->SetUnitCandidates("rad deg");

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorMessenger::~mpdPrimaryGeneratorMessenger()
{
 // delete fRndmCmd;
  delete fGunDir;
  delete thetaCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if(command == thetaCmd)
  {
  fAction->SetTheta(thetaCmd->GetNewDoubleValue(newValue));

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

