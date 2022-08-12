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
/// \file mpdCalorimeterSD.hh
/// \brief Definition of the mpdCalorimeterSD class

#ifndef mpdCalorimeterSD_h
#define mpdCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"

#include "mpdCalorHit.hh"

#include <vector>

class G4Step;
//class G4Track; // new
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class mpdCalorimeterSD : public G4VSensitiveDetector
{
  public:
    mpdCalorimeterSD(const G4String& name,
                     const G4String& hitsCollectionName, 
                     G4int nofCells);
    virtual ~mpdCalorimeterSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
    
    G4int oldID;
    G4int ParentID;  
    G4int eventID;
    G4double kineticPionEnergy;
    G4double MX_pion;
    G4double MY_pion;
    G4double MZ_pion;
    G4double M_Pion_total;
    G4double Theta_pion;
    G4double Phi_pion;
    G4double MX_pion_pass;
    G4double MY_pion_pass;
    G4double MZ_pion_pass;
    G4double X_pion;
    G4double Y_pion;
    G4double Z_pion;
   
/*
    G4double kineticPionEnergy =0.;
    G4double MX_pion =-999999999.;
    G4double MY_pion = -999999999.;
    G4double MZ_pion = -999999999.;
    G4double M_Pion_total = -999999999.;
    G4double Theta_pion = -999999999.;
    G4double Phi_pion = -999999999.;
    G4double MX_pion_pass = -999999999.;
    G4double MY_pion_pass = -999999999.;
    G4double MZ_pion_pass = -999999999.; 
    G4double X_pion = -999999999.;
    G4double Y_pion = -999999999.;
    G4double Z_pion = -999999999.;
    */
  private:
    mpdCalorHitsCollection* fHitsCollection;
    G4int  fNofCells;
   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

