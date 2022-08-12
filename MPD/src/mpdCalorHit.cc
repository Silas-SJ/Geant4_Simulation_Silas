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
/// \file mpdCalorHit.cc
/// \brief Implementation of the mpdCalorHit class

#include "mpdCalorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<mpdCalorHit>* mpdCalorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorHit::mpdCalorHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fLayerID(-111),
   fkineticPionEnergy(-999999999.),
   fMX_pion(-999999999.),
   fMY_pion(-999999999.),
   fMZ_pion(-999999999.),
   fTheta_pion(-999999999.),
   fPhi_pion(-999999999.),
   fX_pion(-999999999.),
   fY_pion(-999999999.),
   fZ_pion(-999999999.),
   fPionDecay(0),
   fMuonDecay(0),
   fPionTrackIDContainer(-111),
   fMuonTrackIDContainer(-111),
   fn_scint(0),
   fn_cerenk(0),
   fenergy_cerenk(0.), 
   ftime_cerenk(0.), 
   fenergy_scint(0.), 
   ftime_scint(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorHit::~mpdCalorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorHit::mpdCalorHit(const mpdCalorHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fLayerID     = right.fLayerID;
  fkineticPionEnergy = right.fkineticPionEnergy;
  fMX_pion = right.fMX_pion;
  fMY_pion = right.fMY_pion;
  fMZ_pion = right.fMZ_pion;
  fTheta_pion = right.fTheta_pion;
  fPhi_pion = right.fPhi_pion;
  fX_pion = right.fX_pion;
  fY_pion = right.fY_pion;
  fZ_pion = right.fZ_pion;
  fPionDecay   = right.fPionDecay;
  fMuonDecay   = right.fMuonDecay;
  fPionTrackIDContainer = right.fPionTrackIDContainer;
  fMuonTrackIDContainer = right.fMuonTrackIDContainer;
  fn_scint   = right.fn_scint;
  fn_cerenk   = right.fn_cerenk;
  fenergy_cerenk = right.fenergy_cerenk;
  ftime_cerenk = right.ftime_cerenk;
  fenergy_scint = right.fenergy_scint;
  ftime_scint = right.ftime_scint;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const mpdCalorHit& mpdCalorHit::operator=(const mpdCalorHit& right)
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fLayerID     = right.fLayerID;
  fkineticPionEnergy = right.fkineticPionEnergy;
  fMX_pion = right.fMX_pion;
  fMY_pion = right.fMY_pion;
  fMZ_pion = right.fMZ_pion;
  fTheta_pion = right.fTheta_pion;
  fPhi_pion = right.fPhi_pion; 
  fX_pion = right.fX_pion;
  fY_pion = right.fY_pion;
  fZ_pion = right.fZ_pion;
  fPionDecay   = right.fPionDecay;
  fMuonDecay   = right.fMuonDecay;
  fPionTrackIDContainer = right.fPionTrackIDContainer;
  fMuonTrackIDContainer = right.fMuonTrackIDContainer;
  fn_scint   = right.fn_scint;
  fn_cerenk   = right.fn_cerenk;
  fenergy_cerenk = right.fenergy_cerenk;
  ftime_cerenk = right.ftime_cerenk;
  fenergy_scint = right.fenergy_scint;
  ftime_scint = right.ftime_scint;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool mpdCalorHit::operator==(const mpdCalorHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdCalorHit::Print()
{
  G4cout
     << "Edep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......