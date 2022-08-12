//
//  TreeManager.cc
//  
//
//  Created by Helio Nogima on 7/11/21.
//

#include "mpdTreeManager.hh"
#include "G4UnitsTable.hh"
#include <TTree.h>
#include <TFile.h>

mpdTreeManager::mpdTreeManager()
:fRootFile(0),
 fTree1(0)
// fEabs(0.), fEgap(0.) ,fLabs(0.), fLgap(0.)
{
  fTree1 = 0;
}

mpdTreeManager::~mpdTreeManager()
{
  if (fRootFile) delete fRootFile;
}

void mpdTreeManager::Book()
{
    G4String fileName = "mpdTree.root";
    fRootFile = new TFile(fileName,"RECREATE");
  //  fRootFile = new TFile(fileName,"UPDATE");
    /*
    if (fRootFile->Get("Ntuple1")) {
       std::cout<< "passa 3" <<  std::endl;
       fTree1 = (TTree*) fRootFile->Get("Ntuple1");
    }else{
      G4cout << " mpdTreeManager::Book :"
             << " creating the ROOT TFile "
             << G4endl;
             std::cout<< "passa 5" <<  std::endl;
     fTree1 = new TTree("Ntuple1", "MPD data");
    }
    */
    fTree1 = new TTree("Ntuple1", "MPD data");
    fTree1->Branch("Theta", &fTheta, "Theta/D");
    fTree1->Branch("Phi", &fPhi, "Phi/D");
    fTree1->Branch("Energy",&fEnergy,"Energy/D");
    fTree1->Branch("Px",&fPx,"Px/D");
    fTree1->Branch("Py",&fPy,"Py/D");
    fTree1->Branch("Pz",&fPz,"Pz/D");
    fTree1->Branch("X",&fX,"X/D");
    fTree1->Branch("Y",&fY,"Y/D");
    fTree1->Branch("Z",&fZ,"Z/D");
                          // For layers /// (Abs)
    fTree1->Branch("AbsPionDecay","std::map<int,int>",&fAbsPionDecay);
    fTree1->Branch("AbsPionPassed","std::map<int,int>",&fAbsPionPassou);
    fTree1->Branch("AbsXPionPosition","std::map<double,double>",&fAbsXPionPosition);
    fTree1->Branch("AbsYPionPosition","std::map<double,double>",&fAbsYPionPosition);
    fTree1->Branch("AbsZPionPosition","std::map<double,double>",&fAbsZPionPosition);
                      
                          //For Scintillator //
    fTree1->Branch("ScintkineticPionEnergy","std::map<double,double>",&fScintkineticPionEnergy);
    fTree1->Branch("ScintMX_Pion","std::map<double,double>",&fScintMX_Pion);
    fTree1->Branch("ScintMY_Pion","std::map<double,double>",&fScintMY_Pion);
    fTree1->Branch("ScintMZ_Pion","std::map<double,double>",&fScintMZ_Pion);
    fTree1->Branch("ScintTheta_Pion","std::map<double,double>",&fScintTheta_Pion); 
    fTree1->Branch("ScintPhi_Pion","std::map<double,double>",&fScintPhi_Pion);
    fTree1->Branch("ScintPionDecay","std::map<int,int>",&fScintPionDecay);
    fTree1->Branch("ScintPionPassed","std::map<int,int>",&fScintPionPassou);
    fTree1->Branch("ScintMuonDecay","std::map<int,int>",&fScintMuonDecay);
    fTree1->Branch("ScintMuonPassed","std::map<int,int>",&fScintMuonPassou);
    fTree1->Branch("ScintNumPhotonScintillation","std::map<int,int>",&fScintn_scint);
    fTree1->Branch("ScintNumPhotonCerenk","std::map<int,int>",&fScintn_cerenk);
    fTree1->Branch("ScintPhotonEnergyCerenk","std::map<double,double>",&fScintenergy_cerenk);
    fTree1->Branch("ScintPhotonTimeCerenk","std::map<double,double>",&fScinttime_cerenk);
    fTree1->Branch("ScintPhotonEnergyScintillation","std::map<double,double>",&fScintenergy_scint);
    fTree1->Branch("ScintPhotonTimeScintillation","std::map<double,double>",&fScinttime_scint);
   
                             //For Detabs //
    fTree1->Branch("DetabsPionDecay","std::map<int,int>",&fDetabsPionDecay);
    fTree1->Branch("DetabsPionPassed","std::map<int,int>",&fDetabsPionPassou);
                               //For Gap //
    fTree1->Branch("GapPionDecay","std::map<int,int>",&fGapPionDecay);
    fTree1->Branch("GapPionPassed","std::map<int,int>",&fGapPionPassou);
    
    fTree1->Branch("OTpmt","std::vector<int>",&fOTpmt);
    fTree1->Branch("OTtime","std::vector<double>",&fOTtime);
//    fTree1->Branch("OTV","std::vector<pair<G4int,G4double>>",&fOTV);
}

void mpdTreeManager::Save()
{
  if (! fRootFile) return;
  
  fRootFile->Write();       // Writing to the file
  fRootFile->Close();       // and closing the tree (and the file)
  
  G4cout << "\n----> Tree saved\n" << G4endl;
}




void mpdTreeManager::FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,G4int> AbsPionDeca, std::map<G4int,G4int> AbsPionPasso, std::map<G4double,G4double> AbsXPionPosi, std::map<G4double,G4double> AbsYPionPosi, std::map<G4double,G4double> AbsZPionPosi, std::map<G4double,G4double> ScintkineticPionEner, std::map<G4double,G4double> ScintMX_pion, std::map<G4double,G4double> ScintMY_pion, std::map<G4double,G4double> ScintMZ_pion, std::map<G4double,G4double> ScintTheta_pion, std::map<G4double,G4double> ScintPhi_pion, std::map<G4int,G4int> ScintPionDeca, std::map<G4int,G4int> ScintPionPasso, std::map<G4int,G4int> ScintMuonDeca, std::map<G4int,G4int> ScintMuonPasso, std::map<G4int,G4int> Scintn_scin, std::map<G4int,G4int> Scintn_ceren, std::map<G4double,G4double> Scintenergy_ceren, std::map<G4double,G4double> Scinttime_ceren, std::map<G4double,G4double> Scintenergy_scin, std::map<G4double,G4double> Scinttime_scin, std::map<G4int,G4int> DetabsPionDeca, std::map<G4int,G4int> DetabsPionPasso, std::map<G4int,G4int> GapPionDeca, std::map<G4int,G4int> GapPionPasso, std::vector<G4int> otp, std::vector<G4double> ott)


{
fTheta = theta;
fPhi = phi;
fEnergy = energy;
fPx = px;
fPy = py;
fPz = pz;
fX = x;
fY = y;
fZ = z;

fOTpmt=otp;
fOTtime=ott;

// For Layers (abs) //
fAbsPionDecay = AbsPionDeca;
fAbsPionPassou = AbsPionPasso;
fAbsXPionPosition = AbsXPionPosi;
fAbsYPionPosition = AbsYPionPosi;
fAbsZPionPosition = AbsZPionPosi;

     // For Scintillator //
fScintkineticPionEnergy = ScintkineticPionEner;
fScintMX_Pion = ScintMX_pion;
fScintMY_Pion = ScintMY_pion;
fScintMZ_Pion = ScintMZ_pion;
fScintTheta_Pion = ScintTheta_pion;
fScintPhi_Pion = ScintPhi_pion;
fScintPionDecay = ScintPionDeca;
fScintPionPassou = ScintPionPasso;
fScintMuonDecay = ScintMuonDeca;
fScintMuonPassou = ScintMuonPasso;
fScintn_scint = Scintn_scin;
fScintn_cerenk = Scintn_ceren;
fScintenergy_cerenk = Scintenergy_ceren;
fScinttime_cerenk = Scinttime_ceren;
fScintenergy_scint = Scintenergy_scin;
fScinttime_scint = Scinttime_scin;
      // For Detabs //
fDetabsPionDecay = DetabsPionDeca;
fDetabsPionPassou = DetabsPionPasso;
   // For Gap //
fGapPionDecay = GapPionDeca;
fGapPionPassou = GapPionPasso;

  if (fTree1) fTree1->Fill();
}
