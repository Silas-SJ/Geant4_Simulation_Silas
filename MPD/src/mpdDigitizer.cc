//
//  mpdDigitizer.cc
//
//
//  Created by Helio Nogima on 7/11/21.
//
#include <TH1D.h>
#include "TSpectrum.h"
#include "TFile.h"
#include "mpdDigitizer.hh"
#include "mpdDigi.hh"
#include "mpdPMTHit.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include <vector>
#include <stdio.h>

mpdDigitizer::mpdDigitizer(G4String name)
  :G4VDigitizerModule(name)
{
  G4String colName = "mpdDigitsCollection";
  collectionName.push_back(colName);
//  G4cout << "Passa no construtor do mpdDigitizer" << G4endl;
}

mpdDigitizer::~mpdDigitizer()
{}

void mpdDigitizer::Digitize()
{
  DigitsCollection = new mpdDigitsCollection
//  	("mpdDigitizer",collectionName[0]); // Criando um Digi Collection
    ("mpdDigitizer","mpdDigitsCollection");
  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();
  G4int PHCID; // pmtHitCollection
  G4int PSHCID; //pmtSignalHitCollection

  
  // PMT Hits collection
  
  PHCID = DigiMan->GetHitsCollectionID("pmtHitCollection");
  mpdPMTHitsCollection* PHC = 0;
  PHC = (mpdPMTHitsCollection*)
    (DigiMan->GetHitsCollection(PHCID));

    PSHCID = DigiMan->GetHitsCollectionID("pmtSignalHitCollection");
  mpdPMTHitsCollection* PSHC = 0;
  PSHC = (mpdPMTHitsCollection*)
     (DigiMan->GetHitsCollection(PSHCID));
/*
G4double Gauss_centro = 13.58*125e-15;
G4double Gauss_sigma = 2.74*125e-15;
G4double Gauss_amplitude = 3.5;


 for(int i=0;i<200;i++){
      probList[i]=(Gauss_amplitude*exp((-1*(pow((i-Gauss_centro),2)))/(pow(Gauss_sigma,2))));
  }
*/

//-----------------------------------------------------------------------------------------------
  if (PHC)
    {
      G4int n_hit = PHC->entries();

      for (G4int i=0;i<n_hit;i++)
        {
	      G4int pmt = (*PHC)[i]->GetPMTNumber();
          G4double energy = (*PHC)[i]->GetEdep();
          G4double PhotonAbs = (*PHC)[i]->GetPhotonAbs();
          if (PhotonAbs){
              G4double Photon_eletron =
                 CLHEP::RandGauss::shoot(13.58,2.74);
              mpdDigi* Digi = new mpdDigi();
              Digi->SetOpticalPhotonEnergy(energy);
              Digi->SetPhoton_eletron(Photon_eletron);
              Digi->SetPMT(pmt);
              Digi->SetOpticalPhotonAbs(PhotonAbs);
         }
        }
    }

//-----------------------------------------------------------------------------------------------
    std::vector <G4double> tsort0;
    std::vector <G4double> tsort1;
    std::vector <G4double> timevec0;
    std::vector <G4double> timevec1;
    std::vector <G4int> pmtvec;
    std::vector <G4double> timevec;
    G4double minTime[]={0,0,0,0,0,0,0,0,0,0,0,0};
    G4double maxTime[]={0,0,0,0,0,0,0,0,0,0,0,0};
    std::map<G4int,std::vector<G4double>> tsort;
    std::map<G4int,std::vector<G4double>> timevector;
    std::map<G4int,G4double> hbegin;
    std::map<G4int,G4double> hend;
    std::map<G4int,G4double> htinterval;
    std::map<G4int,G4int> nbin;
    std::map<G4int,G4double> binwidth;
    G4double hbegin0=0.;
    G4double hend0=0.;
    G4double hbegin1=0.;
    G4double hend1=0.;
    G4int npeaks0=0;
    G4int npeaks1=0; 
   if (PSHC)
        {
          G4int n_hit = PSHC->entries();           
          for (G4int i=0;i<n_hit;i++)
            {

              G4int pmt = (*PSHC)[i]->GetPMTNumber();
              G4double stime = (*PSHC)[i]->GetPMTTime();
              if(pmt==0) tsort0.push_back(stime);
              if(pmt==1) tsort1.push_back(stime);
              tsort[pmt].push_back(stime);
              G4double scharge = (*PSHC)[i]->GetPMTCharge();
                if (i==0) {
                    minTime[pmt]=stime;
                    maxTime[pmt]=stime;
                }else{
                    if (minTime[pmt]>stime) minTime[pmt]=stime;
                    if (maxTime[pmt]<stime) maxTime[pmt]=stime;
                }

            }
            sort(tsort0.begin(),tsort0.end());
            if (tsort0.size()>0) {
                hbegin0=tsort0[0];
                hend0=tsort0[tsort0.size()-1];
            }
            sort(tsort1.begin(),tsort1.end());
            if (tsort1.size()>0) {
                hbegin1=tsort1[0];
                hend1=tsort1[tsort1.size()-1];
	    }
         TFile *fp = new TFile("signal.root","recreate");
         std::map<G4int,TH1F*> htime;
         for (auto i=tsort.begin();i!=tsort.end();++i){
          sort(tsort[i->first].begin(),tsort[i->first].end());
          if (tsort[i->first].size()>0) {
            G4double exTime=10.;
            hbegin[i->first]=tsort[i->first][0]-exTime/2;
            hend[i->first]=tsort[i->first][tsort[i->first].size()-1]+exTime/2;
            htinterval[i->first]=hend[i->first]-hbegin[i->first];
            nbin[i->first]=round(htinterval[i->first]);
            char hname[8];
            sprintf(hname,"h2time%d",i->first);
            htime[i->first] = new TH1F(hname,"",nbin[i->first],hbegin[i->first],hend[i->first]);
            binwidth[i->first]=htinterval[i->first]*1.e-9/htime[i->first]->GetNbinsX();
          }
         }
         G4int nbin0=round((hend0-hbegin0+10));
	 G4int nbin1=round((hend1-hbegin1+10));
         auto htime0 = new TH1F("htime0","PMT0 time",round((hend0-hbegin0+10)),hbegin0-5,hend0+5);
         auto htime1 = new TH1F("htime1","PMT1 time",round((hend1-hbegin1+10)),hbegin1-5,hend1+5);
         G4double interval0 = hend0-hbegin0+10;
         G4double interval1 = hend1-hbegin1+10;
	 G4double bintime0 = interval0*1.e-9/htime0->GetNbinsX(); 
         G4double bintime1 = interval1*1.e-9/htime1->GetNbinsX(); 
         for (G4int i=0;i<n_hit;i++)
           {
             if((*PSHC)[i]->GetEdep()>1.*CLHEP::eV){
              G4int pmt = (*PSHC)[i]->GetPMTNumber();
              G4double stime = (*PSHC)[i]->GetPMTTime();
              G4double scharge = (*PSHC)[i]->GetPMTCharge();
              if (pmt==0) htime0->Fill(stime,50*scharge/bintime0); 
              if (pmt==1) htime1->Fill(stime,50*scharge/bintime1);
              htime[pmt]->Fill(stime,50*scharge/binwidth[pmt]);
             } 
           }
          
           for(int i = 1; i <= htime0->GetNbinsX(); i++ ) if(htime0->GetBinContent(i)>0.015) timevec0.push_back(htime0->GetBinCenter(i));
           for(int i = 1; i <= htime0->GetNbinsX(); i++ ) if(htime1->GetBinContent(i)>0.015) timevec1.push_back(htime1->GetBinCenter(i));

           for(auto j=htime.begin();j!=htime.end();++j){

            for(int i = 1; i <= htime[j->first]->GetNbinsX(); i++ ) if(htime[j->first]->GetBinContent(i)>0.015) timevector[j->first].push_back(htime[j->first]->GetBinCenter(i));
            for(int i = 1; i <= htime[j->first]->GetNbinsX(); i++ ){
              if(htime[j->first]->GetBinContent(i)>0.015){
               pmtvec.push_back(j->first);
               timevec.push_back(htime[j->first]->GetBinCenter(i));
              }
            }
           }
           mpdDigi* SignalDigi = new mpdDigi();
 //          SignalDigi->SetOverThresholdTimeVector0(timevec0);
 //          SignalDigi->SetOverThresholdTimeVector1(timevec1);
 //          SignalDigi->SetOverThresholdTimeVector(timevector);
           SignalDigi->SetOverThresholdPMTVec(pmtvec);
           SignalDigi->SetOverThresholdTimeVec(timevec);
           DigitsCollection->insert(SignalDigi);
           fp->Write("htime0");
           fp->Write("htime1");
           fp->Write("htime");
           fp->Close();
        }
 
//-------------------------------------------------------------------------------------------------

  StoreDigiCollection(DigitsCollection);
 G4int DCID = -1;
  if(DCID<0)
    { 
      	   DigiMan->List();
    DCID = DigiMan->GetDigiCollectionID("mpdDigitizer/DigitsCollection");
       // G4cout << "Digitizer DCID: " << DCID << G4endl;
    }
  
  
}
