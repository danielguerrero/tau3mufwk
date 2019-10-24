#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TChain.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>   
#include "SelectionHelper.h"
#include "MultiplotterHelper.h"
#include "MuonTriggerHelper.h"
#include "MuonTriggerRates.h"
#include <string>
#include "CfgParser.h"

using namespace std;

void TriggerORRateCalculator(std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons, int trig1, int trig2, int trig3)
{
int trigger=0;
//int trigger1plus=0;
for (uint i=0; i<L1TTMuons.size(); ++i)
{  
 std::vector<std::tuple<TLorentzVector,int,float>> triggered1,triggered2,triggered3;
 if(trig1==1 )  triggered1 = PureTriggerSelector(L1TTMuons.at(i)      ,                         0.3, 0.6, 1.0);  
 if(trig1==2 )  triggered1 = PureTriggerSelector(L1TkMuStubMuons.at(i),                         0.3, 0.6, 1.0);
 if(trig1==3 )  triggered1 = PureTriggerSelector(L1TkMuMuons.at(i),                             0.3, 0.6, 1.0);
 if(trig1==4 )  triggered1 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TTMuons.at(i)  ,      0.3, 0.6, 1.0); 
 if(trig1==5 )  triggered1 = MixedTriggerSelector(L1TkMuMuons.at(i),    L1TTMuons.at(i)  ,      0.3, 0.6, 1.0);
 if(trig1==6 )  triggered1 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TkMuMuons.at(i),      0.3, 0.6, 1.0); 
 if(trig1==7 )  triggered1 = MixedTriggerSelector(L1TkMuMuons.at(i),L1TkMuStubMuons.at(i),      0.3, 0.6, 1.0);
 if(trig1==8 )  triggered1 = MixedME0TriggerSelector(L1TTMuons.at(i),ME0StubMuons.at(i),        0.9, 0.6, 1.0);
 if(trig1==9 )  triggered1 = MixedME0TriggerSelector(L1TkMuStubMuons.at(i),ME0StubMuons.at(i),  0.9, 0.6, 1.0);
 if(trig1==10)  triggered1 = MixedME0TriggerSelector(L1TkMuMuons.at(i),ME0StubMuons.at(i),      0.9, 0.6, 1.0);
 if(trig2==1 )  triggered2 = PureTriggerSelector(L1TTMuons.at(i)      ,                         0.3, 0.6, 1.0);  
 if(trig2==2 )  triggered2 = PureTriggerSelector(L1TkMuStubMuons.at(i),                         0.3, 0.6, 1.0);
 if(trig2==3 )  triggered2 = PureTriggerSelector(L1TkMuMuons.at(i),                             0.3, 0.6, 1.0);
 if(trig2==4 )  triggered2 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TTMuons.at(i)  ,      0.3, 0.6, 1.0); 
 if(trig2==5 )  triggered2 = MixedTriggerSelector(L1TkMuMuons.at(i),    L1TTMuons.at(i)  ,      0.3, 0.6, 1.0);
 if(trig2==6 )  triggered2 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TkMuMuons.at(i),      0.3, 0.6, 1.0); 
 if(trig2==7 )  triggered2 = MixedTriggerSelector(L1TkMuMuons.at(i),L1TkMuStubMuons.at(i),      0.3, 0.6, 1.0);
 if(trig2==8 )  triggered2 = MixedME0TriggerSelector(L1TTMuons.at(i),ME0StubMuons.at(i),        0.9, 0.6, 1.0);
 if(trig2==9 )  triggered2 = MixedME0TriggerSelector(L1TkMuStubMuons.at(i),ME0StubMuons.at(i),  0.9, 0.6, 1.0);
 if(trig2==10)  triggered2 = MixedME0TriggerSelector(L1TkMuMuons.at(i),ME0StubMuons.at(i),      0.9, 0.6, 1.0);
 if(trig3==1 )  triggered3 = PureTriggerSelector(L1TTMuons.at(i)      ,                         0.3, 0.6, 1.0);  
 if(trig3==2 )  triggered3 = PureTriggerSelector(L1TkMuStubMuons.at(i),                         0.3, 0.6, 1.0);
 if(trig3==3 )  triggered3 = PureTriggerSelector(L1TkMuMuons.at(i),                             0.3, 0.6, 1.0);
 if(trig3==4 )  triggered3 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TTMuons.at(i)  ,      0.3, 0.6, 1.0); 
 if(trig3==5 )  triggered3 = MixedTriggerSelector(L1TkMuMuons.at(i),    L1TTMuons.at(i)  ,      0.3, 0.6, 1.0);
 if(trig3==6 )  triggered3 = MixedTriggerSelector(L1TkMuStubMuons.at(i),L1TkMuMuons.at(i),      0.3, 0.6, 1.0); 
 if(trig3==7 )  triggered3 = MixedTriggerSelector(L1TkMuMuons.at(i),L1TkMuStubMuons.at(i),      0.3, 0.6, 1.0);
 if(trig3==8 )  triggered3 = MixedME0TriggerSelector(L1TTMuons.at(i),ME0StubMuons.at(i),        0.9, 0.6, 1.0);
 if(trig3==9 )  triggered3 = MixedME0TriggerSelector(L1TkMuStubMuons.at(i),ME0StubMuons.at(i),  0.9, 0.6, 1.0);
 if(trig3==10)  triggered3 = MixedME0TriggerSelector(L1TkMuMuons.at(i),ME0StubMuons.at(i),      0.9, 0.6, 1.0);
 if(triggered1.size() == 3 || triggered2.size() == 3 || triggered3.size()==3 ) trigger++;
}
//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Acceptance for Triggers %i, %i, %i    ",trig1,trig2,trig3)<<(float)trigger/L1TTMuons.size()<<endl;
cout<<Form("Maximum Rate for Triggers %i, %i, %i  (kHz) = ",trig1,trig2,trig3)<<(float)2760*11.246*trigger/L1TTMuons.size()<<endl;
}


void TriggerRateStudies(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons)
{
//Set three track triggers (cut on charge,                   mass,  dR,  dZ)
PureTripleMuonTriggerCalculator(L1TTMuons                   , 0.3, 0.6, 1.0, "1");
PureTripleMuonTriggerCalculator(L1TkMuStubMuons             , 0.3, 0.6, 1.0, "2");
PureTripleMuonTriggerCalculator(L1TkMuMuons                 , 0.3, 0.6, 1.0, "3");
MixedTripleMuonTriggerCalculator(L1TkMuStubMuons,L1TTMuons  , 0.3, 0.6, 1.0, "4");
MixedTripleMuonTriggerCalculator(L1TkMuMuons,L1TTMuons      , 0.3, 0.6, 1.0, "5");
MixedTripleMuonTriggerCalculator(L1TkMuStubMuons,L1TkMuMuons, 0.3, 0.6, 1.0, "6");
MixedTripleMuonTriggerCalculator(L1TkMuMuons,L1TkMuStubMuons, 0.3, 0.6, 1.0, "7");
//Triggers using two tracks + ME0 stub (cut on charge,        dZ,  dR,  ME0 bending)
MixedTripleME0TriggerCalculator(L1TTMuons,      ME0StubMuons, 0.9, 0.6, 1.0, "8");
MixedTripleME0TriggerCalculator(L1TkMuStubMuons,ME0StubMuons, 0.9, 0.6, 1.0, "9");
MixedTripleME0TriggerCalculator(    L1TkMuMuons,ME0StubMuons, 0.9, 0.6, 1.0, "10");
//Make OR calculator
TriggerORRateCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons, 2, 7,   7);
TriggerORRateCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons, 2, 10, 10);
TriggerORRateCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons, 7, 10, 10);
TriggerORRateCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons, 2,  7, 10);
}  

void MuonTriggerRates(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons, 
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons,
	string outputname)
{
cout << "[INFO] MuonTriggerRates module running ..."<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_muontriggerrates.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting Trigger Studies"<<endl;
TriggerRateStudies(EMTFMuons,L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons);

tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  
