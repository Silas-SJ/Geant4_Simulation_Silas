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
/// \file mpdCalorHit.hh
/// \brief Definition of the mpdCalorHit class

#ifndef mpdCalorHit_h
#define mpdCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include <set> 

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class mpdCalorHit : public G4VHit
{
  public:
    mpdCalorHit();
    mpdCalorHit(const mpdCalorHit&);
    virtual ~mpdCalorHit();

    // operators
    const mpdCalorHit& operator=(const mpdCalorHit&);
    G4bool operator==(const mpdCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    //void Add(G4double de, G4double dl, G4int layer, G4double kineticPionEnergy, G4double MX_pion, G4double MY_pion, G4double MZ_pion, G4double Theta_pion, G4double Phi_pion, G4double X_pion, G4double Y_pion, G4double Z_pion, G4int n_scint, G4int n_cerenk, G4double energy_cerenk, G4double time_cerenk, G4double energy_scint, G4double time_scint);
    void Add(G4double de, G4double dl, G4int layer, G4int n_scint, G4int n_cerenk, G4double energy_cerenk, G4double time_cerenk, G4double energy_scint, G4double time_scint);
    void AddPionDecay();
    void AddMuonDecay();
    void SetPionTrackID(G4int tid);
    void SetMuonTrackID(G4int tid);
    void SetX_PionPosition(G4double tid);
    void SetY_PionPosition(G4double tid);
    void SetZ_PionPosition(G4double tid);
    void SetkineticPionEnergy(G4double tid);
    void SetMXPion(G4double tid);
    void SetMYPion(G4double tid);
    void SetMZPion(G4double tid);
    void SetThetaPion(G4double tid);
    void SetPhiPion(G4double tid);
    // get methods
    G4double GetEdep() ;
    G4double GetTrackLength() ;
    G4int GetLayerID() ;
    G4double GetkineticPionEnergy() ;
    G4double GetMXPion() ;
    G4double GetMYPion() ;
    G4double GetMZPion() ;
    G4double GetThetaPion() ;
    G4double GetPhiPion() ;
    G4double GetXPionPosition() ;
    G4double GetYPionPosition() ;
    G4double GetZPionPosition() ;
    G4int GetPionDecay() ;
    G4int GetMuonDecay() ;  
    G4int GetPionTrackID() ;    
    G4int GetMuonTrackID() ;
    G4int GetNumberScintillation() ; 
    G4int GetNumberCerenkov() ;
    G4double GetEnergyCerenkov() ; 
    G4double GetTimeCerenkov() ;
    G4double GetEnergyScintillation() ;
    G4double GetTimeScintillation() ;

  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    G4int fLayerID;
    G4double fkineticPionEnergy;
    G4double fMX_pion;
    G4double fMY_pion;
    G4double fMZ_pion;
    G4double fTheta_pion;
    G4double fPhi_pion;
    G4double fX_pion;
    G4double fY_pion;
    G4double fZ_pion;
    G4int fPionDecay;
    G4int fMuonDecay;
    G4int fPionTrackIDContainer;
    G4int fMuonTrackIDContainer;
    G4int fn_scint;
    G4int fn_cerenk;
    G4double fenergy_cerenk;
    G4double ftime_cerenk;
    G4double fenergy_scint;
    G4double ftime_scint;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using mpdCalorHitsCollection = G4THitsCollection<mpdCalorHit>;

extern G4ThreadLocal G4Allocator<mpdCalorHit>* mpdCalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* mpdCalorHit::operator new(size_t)
{
  if (!mpdCalorHitAllocator) {
    mpdCalorHitAllocator = new G4Allocator<mpdCalorHit>;
  }
  void *hit;
  hit = (void *) mpdCalorHitAllocator->MallocSingle();
  return hit;
}

inline void mpdCalorHit::operator delete(void *hit)
{
  if (!mpdCalorHitAllocator) {
    mpdCalorHitAllocator = new G4Allocator<mpdCalorHit>;
  }
  mpdCalorHitAllocator->FreeSingle((mpdCalorHit*) hit);
}

//inline void mpdCalorHit::Add(G4double de, G4double dl, G4int layer, G4double kineticPionEnergy,  G4double MX_pion, G4double MY_pion, G4double MZ_pion, G4double Theta_pion, G4double Phi_pion, G4double X_pion, G4double Y_pion, G4double Z_pion, G4int n_scint, G4int n_cerenk, G4double energy_cerenk, G4double time_cerenk, G4double energy_scint, G4double time_scint) {
inline void mpdCalorHit::Add(G4double de, G4double dl, G4int layer, G4int n_scint, G4int n_cerenk, G4double energy_cerenk, G4double time_cerenk, G4double energy_scint, G4double time_scint) {
  fEdep += de;
  fTrackLength += dl;
  fLayerID = layer;
  /*
  fkineticPionEnergy = kineticPionEnergy;
  fMX_pion = MX_pion;
  fMY_pion = MY_pion;
  fMZ_pion = MZ_pion;
  fTheta_pion = Theta_pion;
  fPhi_pion = Phi_pion;
  fX_pion = X_pion;
  fY_pion = Y_pion;
  fZ_pion = Z_pion;
  */
  fn_scint = n_scint;
  fn_cerenk = n_cerenk;
  fenergy_cerenk = energy_cerenk;
  ftime_cerenk = time_cerenk;
  fenergy_scint = energy_scint;
  ftime_scint = time_scint;
}
inline void mpdCalorHit::AddPionDecay() {
  fPionDecay++;
}
inline void mpdCalorHit::AddMuonDecay() {
  fMuonDecay++;
}
inline void mpdCalorHit::SetPionTrackID(G4int tid){
//fPionTrackIDContainer.insert(tid);
fPionTrackIDContainer=tid;
}
inline void mpdCalorHit::SetMuonTrackID(G4int tid){
//fMuonTrackIDContainer.insert(tid);
fMuonTrackIDContainer=tid;
}
inline void mpdCalorHit::SetX_PionPosition(G4double tid){
fX_pion=tid;
}
inline void mpdCalorHit::SetY_PionPosition(G4double tid){
fY_pion=tid;
}
inline void mpdCalorHit::SetZ_PionPosition(G4double tid){
fZ_pion=tid;
}
inline void mpdCalorHit::SetkineticPionEnergy(G4double tid){
fkineticPionEnergy=tid;
}
inline void mpdCalorHit::SetMXPion(G4double tid){
fMX_pion=tid;
}
inline void mpdCalorHit::SetMYPion(G4double tid){
fMY_pion=tid;
}
inline void mpdCalorHit::SetMZPion(G4double tid){
fMZ_pion=tid;
}
inline void mpdCalorHit::SetThetaPion(G4double tid){
fTheta_pion=tid;
}
inline void mpdCalorHit::SetPhiPion(G4double tid){
fPhi_pion=tid;
}

inline G4double mpdCalorHit::GetEdep() {
  return fEdep;
}
inline G4double mpdCalorHit::GetTrackLength() {
  return fTrackLength;
}
inline G4int mpdCalorHit::GetLayerID(){
 return fLayerID;
}
inline G4double mpdCalorHit::GetkineticPionEnergy()  {
  return fkineticPionEnergy;
}
inline G4double mpdCalorHit::GetMXPion()  {
  return fMX_pion;
}
inline G4double mpdCalorHit::GetMYPion() {
  return fMY_pion;
}
inline G4double mpdCalorHit::GetMZPion()  {
  return fMZ_pion;
}
inline G4double mpdCalorHit::GetThetaPion()  {
  return fTheta_pion;
}
inline G4double mpdCalorHit::GetPhiPion()  {
  return fPhi_pion;
}
inline G4int mpdCalorHit::GetPionDecay()  {
  return fPionDecay;
}
inline G4int mpdCalorHit::GetMuonDecay()  {
  return fMuonDecay;
}
inline G4int mpdCalorHit::GetPionTrackID()  {
 return fPionTrackIDContainer;
}
inline G4int mpdCalorHit::GetMuonTrackID()  {
 return fMuonTrackIDContainer;
}
inline G4double mpdCalorHit::GetXPionPosition()  {
  return fX_pion;
}
inline G4double mpdCalorHit::GetYPionPosition()  {
  return fY_pion;
}
inline G4double mpdCalorHit::GetZPionPosition()  {
  return fZ_pion;
}
inline G4int mpdCalorHit::GetNumberScintillation()  {
 return fn_scint;
}
inline G4int mpdCalorHit::GetNumberCerenkov()  {
 return fn_cerenk;
}
inline G4double mpdCalorHit::GetEnergyCerenkov()  {
 return fenergy_cerenk;
}
inline G4double mpdCalorHit::GetTimeCerenkov()  {
 return ftime_cerenk;
}
inline G4double mpdCalorHit::GetEnergyScintillation()  {
 return fenergy_scint;
}
inline G4double mpdCalorHit::GetTimeScintillation()  {
 return ftime_scint;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
