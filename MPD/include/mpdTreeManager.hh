//
//  mpdTreeManager.hh
//  
//
//  Created by Helio Nogima on 7/11/21.
//

#ifndef mpdTreeManager_hh
#define mpdTreeManager_hh

#include <stdio.h>
#include "globals.hh"
#include <vector>

class TFile;
class TTree;

class mpdTreeManager
{
  public:
    mpdTreeManager();
    ~mpdTreeManager();

    void Book();
    void Save();
//    void FillNtuple(G4double energyAbs, G4double energyGap, G4double trackLAbs, G4double trackLGap, G4int Npe, G4double PMTCharge1[12], G4double PMTCharge2[12], G4double PMTMediumTime1[12], G4double PMTMediumTime2[12], G4double ScintCharge[12]);
//    void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<G4double> ottv0, std::vector<G4double> ottv1);
//   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,std::vector<G4double>> ottv);  

/*
void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,G4int> AbsPionDecay, std::map<G4int,G4int> AbsPionPassou, std::map<G4int,G4int> ScintPionDecay, std::map<G4int,G4int> ScintPionPassou, std::map<G4int,G4int> ScintMuonDecay, std::map<G4int,G4int> ScintMuonPassou, std::map<G4int,G4int> DetabsPionDecay, std::map<G4int,G4int> DetabsPionPassou, std::map<G4int,G4int> GapPionDecay, std::map<G4int,G4int> GapPionPassou, std::vector<G4int> ovp, std::vector<G4double> ott); 
*/

   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,G4int> AbsPionDecay, std::map<G4int,G4int> AbsPionPassou, std::map<G4double,G4double> AbsXPionPosition, std::map<G4double,G4double> AbsYPionPosition, std::map<G4double,G4double> AbsZPionPosition, std::map<G4double,G4double> ScintkineticPionEnergy, std::map<G4double,G4double> ScintMX_Pion, std::map<G4double,G4double> ScintMY_Pion, std::map<G4double,G4double> ScintMZ_Pion, std::map<G4double,G4double> ScintTheta_Pion, std::map<G4double,G4double> ScintPhi_Pion, std::map<G4int,G4int> ScintPionDecay, std::map<G4int,G4int> ScintPionPassou, std::map<G4int,G4int> ScintMuonDecay, std::map<G4int,G4int> ScintMuonPassou, std::map<G4int,G4int> Scintn_scint, std::map<G4int,G4int> Scintn_cerenk, std::map<G4double,G4double> Scintenergy_cerenk, std::map<G4double,G4double> Scinttime_cerenk, std::map<G4double,G4double> Scintenergy_scint, std::map<G4double,G4double> Scinttime_scint, std::map<G4int,G4int> DetabsPionDecay, std::map<G4int,G4int> DetabsPionPassou, std::map<G4int,G4int> GapPionDecay, std::map<G4int,G4int> GapPionPassou, std::vector<G4int> ovp, std::vector<G4double> ott);  
 
 /*  
 void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, G4double scintPionEnergy, G4double scintPionMX ,G4double scintPionMY, G4double scintPionMZ, std::map<G4int,G4int> AbsPionDecay, std::map<G4int,G4int> AbsPionPassou, std::map<G4int,G4int> ScintPionDecay, std::map<G4int,G4int> ScintPionPassou, std::map<G4int,G4int> ScintMuonDecay, std::map<G4int,G4int> ScintMuonPassou, std::map<G4int,G4int> DetabsPionDecay, std::map<G4int,G4int> DetabsPionPassou, std::map<G4int,G4int> GapPionDecay, std::map<G4int,G4int> GapPionPassou, std::vector<G4int> ovp, std::vector<G4double> ott); 
   */
//   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<std::pair<int,double>> otv);  
  private:
    TFile*   fRootFile;
    TTree*   fTree1;

    G4double fTheta;
    G4double fPhi;
    G4double fEnergy;
    G4double fPx, fPy, fPz;
    G4double fX, fY, fZ;
 /*
    G4double fEabs;
    G4double fEgap;
    G4double fLabs;
    G4double fLgap;
    */
    
   /*
    G4int fPionDecay;
    G4int fMuonDecay;
    G4int fPionPassed;
    G4int fMuonPassed;
  */ 
             // for Layers (Abs) ///
    std::map<G4int,G4int> fAbsPionDecay;
    std::map<G4int,G4int> fAbsPionPassou;
    std::map<G4double,G4double> fAbsXPionPosition;
    std::map<G4double,G4double> fAbsYPionPosition;
    std::map<G4double,G4double> fAbsZPionPosition;
  //  std::map<G4double,G4int> fAbsTheta_Pion_Passed;
   // std::map<G4double,G4int> fAbsPhi_Pion_Passed;
    
              // for Scintillator ///
    std::map<G4double,G4double> fScintkineticPionEnergy;
    std::map<G4double,G4double> fScintMX_Pion;
    std::map<G4double,G4double> fScintMY_Pion;
    std::map<G4double,G4double> fScintMZ_Pion;
    std::map<G4double,G4double> fScintTheta_Pion;
    std::map<G4double,G4double> fScintPhi_Pion;
    
    std::map<G4int,G4int> fScintPionDecay;
    std::map<G4int,G4int> fScintPionPassou;
    std::map<G4int,G4int> fScintMuonDecay;
    std::map<G4int,G4int> fScintMuonPassou;
    std::map<G4int,G4int> fScintn_scint;
    std::map<G4int,G4int> fScintn_cerenk;
    std::map<G4double,G4double> fScintenergy_cerenk;
    std::map<G4double,G4double> fScinttime_cerenk; 
    std::map<G4double,G4double> fScintenergy_scint; 
    std::map<G4double,G4double> fScinttime_scint;
            // for Detabs ///
    std::map<G4int,G4int> fDetabsPionDecay;
    std::map<G4int,G4int> fDetabsPionPassou;  
             // for Gap ///
    std::map<G4int,G4int> fGapPionDecay;
    std::map<G4int,G4int> fGapPionPassou;
    
    std::vector <G4int> fOTpmt;
    std::vector <G4double> fOTtime;
//    std::vector <std::pair<G4int,G4double>> fOTV;
//    std::vector <G4double> fOTTimeVector0;
//    std::vector <G4double> fOTTimeVector1;
//    std::map<G4int,std::vector<G4double>> fOTTimeVector;
  
 //  G4int    fNpe;
//    G4double fPMTCharge1[12], fPMTCharge2[12];
//    G4double fPMTMediumTime1[12], fPMTMediumTime2[12];
//    G4double fScintCharge[12];
};

#endif /* mpdTreeManager_hh */
