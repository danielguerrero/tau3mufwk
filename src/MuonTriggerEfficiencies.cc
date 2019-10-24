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
#include "MuonTriggerHelper.h"
#include "MuonTriggerEfficiencies.h"
#include "MultiplotterHelper.h"
#include <string>
#include "CfgParser.h"

using namespace std;

void TriggerOREfficiencyCalculator(std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
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
cout<<Form("Trigger Efficiency for Triggers %i, %i, %i = ",trig1,trig2,trig3)<<(float)trigger/L1TTMuons.size()<<endl;
}

void TriggerEfficiencyStudies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> EMTFMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons)
{
//Book Histogram for the reference
//MAKE BINS
Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0};
int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
//Make Histograms for the reference
TH1F *h_REF_gen_mu_p_0     =  new TH1F("h_REF_gen_mu_p_0", "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_0    =  new TH1F("h_REF_gen_mu_pt_0", "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_0   =  new TH1F("h_REF_gen_mu_eta_0","", netabins, etabins);
TH1F *h_REF_gen_mu_p_1     =  new TH1F("h_REF_gen_mu_p_1", "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_1    =  new TH1F("h_REF_gen_mu_pt_1", "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_1   =  new TH1F("h_REF_gen_mu_eta_1","", netabins, etabins);
TH1F *h_REF_gen_mu_p_2     =  new TH1F("h_REF_gen_mu_p_2", "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_2    =  new TH1F("h_REF_gen_mu_pt_2", "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_2   =  new TH1F("h_REF_gen_mu_eta_2","", netabins, etabins);
TH1F *h_REF_gen_mu_p_min   =  new TH1F("h_REF_gen_mu_p_min", "", npbins,pbins);
TH1F *h_REF_gen_mu_eta_max =  new TH1F("h_REF_gen_mu_eta_max", "", netabins,etabins);

for (uint i=0; i<GenMuons.size(); ++i)
{  
  h_REF_gen_mu_p_0->Fill(   GenMuons.at(i).at(0).P());
  h_REF_gen_mu_pt_0->Fill(  GenMuons.at(i).at(0).Pt());
  h_REF_gen_mu_eta_0->Fill( abs(GenMuons.at(i).at(0).Eta()) );
  h_REF_gen_mu_p_1->Fill(   GenMuons.at(i).at(1).P());
  h_REF_gen_mu_pt_1->Fill(  GenMuons.at(i).at(1).Pt());
  h_REF_gen_mu_eta_1->Fill( abs(GenMuons.at(i).at(1).Eta()) );
  h_REF_gen_mu_p_2->Fill(   GenMuons.at(i).at(2).P());
  h_REF_gen_mu_pt_2->Fill(  GenMuons.at(i).at(2).Pt());
  h_REF_gen_mu_eta_2->Fill( abs(GenMuons.at(i).at(2).Eta()) );
  h_REF_gen_mu_p_min->Fill( std::min( { GenMuons.at(i).at(0).P(), GenMuons.at(i).at(1).P(), GenMuons.at(i).at(2).P() }));
  h_REF_gen_mu_eta_max->Fill( std::max( { abs(GenMuons.at(i).at(0).Eta()), abs(GenMuons.at(i).at(1).Eta()), abs(GenMuons.at(i).at(2).Eta()) })); 
}

//Check matching between tracks and genmuons
MatchingCalculator(GenMuons,L1TTMuons,"Tk");
MatchingCalculator(GenMuons,L1TkMuStubMuons,"TkMuStub");
MatchingCalculator(GenMuons,L1TkMuMuons,"TkMu");

//Set three track triggers (cut on charge,                         mass, dR,  dZ)
PureMuonTriggerCalculator(GenMuons,L1TTMuons,                      0.3, 0.6, 1.0, "1");
PureMuonTriggerCalculator(GenMuons,L1TkMuStubMuons,                0.3, 0.6, 1.0, "2");
PureMuonTriggerCalculator(GenMuons,L1TkMuMuons,                    0.3, 0.6, 1.0, "3");
MixedMuonTriggerCalculator(GenMuons,L1TkMuStubMuons,L1TTMuons,     0.3, 0.6, 1.0, "4");
MixedMuonTriggerCalculator(GenMuons,L1TkMuMuons,L1TTMuons,         0.3, 0.6, 1.0, "5");
MixedMuonTriggerCalculator(GenMuons,L1TkMuStubMuons,L1TkMuMuons,   0.3, 0.6, 1.0, "6");
MixedMuonTriggerCalculator(GenMuons,L1TkMuMuons,L1TkMuStubMuons,   0.3, 0.6, 1.0, "7");
//Triggers using two tracks + ME0 stub (cut on charge,              dz,  dR, ME0 bending)
MixedME0TriggerCalculator(GenMuons,L1TTMuons,       ME0StubMuons,  0.9, 0.6, 1.0, "8");
MixedME0TriggerCalculator(GenMuons,L1TkMuStubMuons, ME0StubMuons,  0.9, 0.6, 1.0, "9");
MixedME0TriggerCalculator(GenMuons,L1TkMuMuons    , ME0StubMuons,  0.9, 0.6, 1.0, "10");
//Make OR calculator
TriggerOREfficiencyCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons,2, 7, 7);
TriggerOREfficiencyCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons,2,10,10);
TriggerOREfficiencyCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons,7,10,10);
TriggerOREfficiencyCalculator(L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons,2, 7,10);
}  

void MixedMuonTriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float mwindow, float drwindow, float dzwindow, string tagname)
{
//Book Histogram for the reference
//MAKE BINS
Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0};
int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
//Make Histograms for the reference
TH1F *h_TRG1_gen_mu_p_0     =  new TH1F(Form("h_TRG%s_gen_mu_p_0",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG1_gen_mu_pt_0    =  new TH1F(Form("h_TRG%s_gen_mu_pt_0",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG1_gen_mu_eta_0   =  new TH1F(Form("h_TRG%s_gen_mu_eta_0",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG1_gen_mu_p_1     =  new TH1F(Form("h_TRG%s_gen_mu_p_1",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG1_gen_mu_pt_1    =  new TH1F(Form("h_TRG%s_gen_mu_pt_1",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG1_gen_mu_eta_1   =  new TH1F(Form("h_TRG%s_gen_mu_eta_1",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG1_gen_mu_p_2     =  new TH1F(Form("h_TRG%s_gen_mu_p_2",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG1_gen_mu_pt_2    =  new TH1F(Form("h_TRG%s_gen_mu_pt_2",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG1_gen_mu_eta_2   =  new TH1F(Form("h_TRG%s_gen_mu_eta_2",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG1_gen_mu_p_min   =  new TH1F(Form("h_TRG%s_gen_mu_p_min",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG1_gen_mu_eta_max =  new TH1F(Form("h_TRG%s_gen_mu_eta_max",tagname.c_str() ), "", netabins,etabins);
TH1F *h_TRG1_invmass        =  new TH1F(Form("h_TRG%s_invmass",tagname.c_str() ), "",100, 0, 3);
TH1F *h_TRG1_mu_p_0         =  new TH1F(Form("h_TRG%s_mu_p_0",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG1_mu_pt_0        =  new TH1F(Form("h_TRG%s_mu_pt_0",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG1_mu_eta_0       =  new TH1F(Form("h_TRG%s_mu_eta_0",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG1_mu_p_1         =  new TH1F(Form("h_TRG%s_mu_p_1",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG1_mu_pt_1        =  new TH1F(Form("h_TRG%s_mu_pt_1",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG1_mu_eta_1       =  new TH1F(Form("h_TRG%s_mu_eta_1",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG1_mu_p_2         =  new TH1F(Form("h_TRG%s_mu_p_2",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG1_mu_pt_2        =  new TH1F(Form("h_TRG%s_mu_pt_2",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG1_mu_eta_2       =  new TH1F(Form("h_TRG%s_mu_eta_2",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG1_variable       =  new TH1F(Form("h_TRG%s_variable",tagname.c_str() ), "", 304, 2, 40);
TH1F *h_TRG1_maxdR          =  new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "", 100, 0,  4);
TH1F *h_TRG1_totalcharge    =  new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG1_mu_dz0_01      =  new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG1_mu_dz0_12      =  new TH1F(Form("h_TRG%s_mu_dz0_12",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG1_mu_dz0_02      =  new TH1F(Form("h_TRG%s_mu_dz0_02",tagname.c_str() ), "", 80, 0, 20);
TH2F *h_TRG1_respt          =  new TH2F(Form("h_TRG%s_respt"    ,tagname.c_str() ), "", 50, 0, 50, 50, 0, 50);

int trigger1=0;
//int trigger1plus=0;
for (uint i=0; i<GenMuons.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = MixedTriggerSelector(Tracks12.at(i),Tracks3.at(i), mwindow, drwindow, dzwindow); 
  if(triggered.size() == 3){
     trigger1++;
     h_TRG1_gen_mu_p_0->Fill(   GenMuons.at(i).at(0).P());
     h_TRG1_gen_mu_pt_0->Fill(  GenMuons.at(i).at(0).Pt());
     h_TRG1_gen_mu_eta_0->Fill( abs(GenMuons.at(i).at(0).Eta()) );
     h_TRG1_gen_mu_p_1->Fill(   GenMuons.at(i).at(1).P());
     h_TRG1_gen_mu_pt_1->Fill(  GenMuons.at(i).at(1).Pt());
     h_TRG1_gen_mu_eta_1->Fill( abs(GenMuons.at(i).at(1).Eta()) );
     h_TRG1_gen_mu_p_2->Fill(   GenMuons.at(i).at(2).P());
     h_TRG1_gen_mu_pt_2->Fill(  GenMuons.at(i).at(2).Pt());
     h_TRG1_gen_mu_eta_2->Fill( abs(GenMuons.at(i).at(2).Eta()) );
     h_TRG1_gen_mu_p_min->Fill( std::min( { GenMuons.at(i).at(0).P(), GenMuons.at(i).at(1).P(), GenMuons.at(i).at(2).P() }));
     h_TRG1_gen_mu_eta_max->Fill( std::max( { abs(GenMuons.at(i).at(0).Eta()), abs(GenMuons.at(i).at(1).Eta()), abs(GenMuons.at(i).at(2).Eta()) })); 
     h_TRG1_invmass->Fill( (get<0>(triggered.at(0)) + get<0>(triggered.at(1)) + get<0>(triggered.at(2)) ).M() );
     h_TRG1_mu_p_0->Fill(   get<0>(triggered.at(0)).P());
     h_TRG1_mu_pt_0->Fill(  get<0>(triggered.at(0)).Pt());
     h_TRG1_mu_eta_0->Fill( abs(get<0>(triggered.at(0)).Eta()) );
     h_TRG1_mu_p_1->Fill(   get<0>(triggered.at(1)).P());
     h_TRG1_mu_pt_1->Fill(  get<0>(triggered.at(1)).Pt());
     h_TRG1_mu_eta_1->Fill( abs(get<0>(triggered.at(1)).Eta()) );
     h_TRG1_mu_p_2->Fill(   get<0>(triggered.at(2)).P());
     h_TRG1_mu_pt_2->Fill(  get<0>(triggered.at(2)).Pt());
     h_TRG1_mu_eta_2->Fill( abs(get<0>(triggered.at(2)).Eta())  );  
     h_TRG1_variable->Fill(  get<0>(triggered.at(1)).Pt());
     h_TRG1_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
     h_TRG1_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
     h_TRG1_respt->Fill( (GenMuons.at(i).at(0)+GenMuons.at(i).at(1)+GenMuons.at(i).at(2)).Pt(), (get<0>(triggered.at(0)) + get<0>(triggered.at(1)) + get<0>(triggered.at(2)) ).Pt() );
     h_TRG1_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
     h_TRG1_mu_dz0_12->Fill( abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) )  );
     h_TRG1_mu_dz0_02->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )  ); 
  }
  else
  {
      h_TRG1_variable->Fill(-1);    
  }

}
MakeRatePlot(h_TRG1_variable,1,"roc");

//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Efficiency for Trigger %s using %.1f GeV window =",tagname.c_str(),mwindow)<<(float)trigger1/GenMuons.size()<<endl;
}

void PureMuonTriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks, float mwindow, float drwindow, float dzwindow, string tagname)
{
//Book Histogram for the reference
//MAKE BINS
Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0};
int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
//Make Histograms for the reference
TH1F *h_TRG2_gen_mu_p_0     =  new TH1F(Form("h_TRG%s_gen_mu_p_0",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG2_gen_mu_pt_0    =  new TH1F(Form("h_TRG%s_gen_mu_pt_0",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG2_gen_mu_eta_0   =  new TH1F(Form("h_TRG%s_gen_mu_eta_0",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG2_gen_mu_p_1     =  new TH1F(Form("h_TRG%s_gen_mu_p_1",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG2_gen_mu_pt_1    =  new TH1F(Form("h_TRG%s_gen_mu_pt_1",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG2_gen_mu_eta_1   =  new TH1F(Form("h_TRG%s_gen_mu_eta_1",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG2_gen_mu_p_2     =  new TH1F(Form("h_TRG%s_gen_mu_p_2",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG2_gen_mu_pt_2    =  new TH1F(Form("h_TRG%s_gen_mu_pt_2",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG2_gen_mu_eta_2   =  new TH1F(Form("h_TRG%s_gen_mu_eta_2",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG2_gen_mu_p_min   =  new TH1F(Form("h_TRG%s_gen_mu_p_min",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG2_gen_mu_eta_max =  new TH1F(Form("h_TRG%s_gen_mu_eta_max",tagname.c_str() ), "", netabins,etabins);
TH1F *h_TRG2_invmass        =  new TH1F(Form("h_TRG%s_invmass",tagname.c_str() ), "",100, 0, 3);
TH1F *h_TRG2_mu_p_0         =  new TH1F(Form("h_TRG%s_mu_p_0",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG2_mu_pt_0        =  new TH1F(Form("h_TRG%s_mu_pt_0",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG2_mu_eta_0       =  new TH1F(Form("h_TRG%s_mu_eta_0",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG2_mu_p_1         =  new TH1F(Form("h_TRG%s_mu_p_1",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG2_mu_pt_1        =  new TH1F(Form("h_TRG%s_mu_pt_1",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG2_mu_eta_1       =  new TH1F(Form("h_TRG%s_mu_eta_1",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG2_mu_p_2         =  new TH1F(Form("h_TRG%s_mu_p_2",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG2_mu_pt_2        =  new TH1F(Form("h_TRG%s_mu_pt_2",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG2_mu_eta_2       =  new TH1F(Form("h_TRG%s_mu_eta_2",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG2_variable       =  new TH1F(Form("h_TRG%s_variable",tagname.c_str() ), "",304, 2, 40);
TH1F *h_TRG2_maxdR          =  new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "", 100, 0,  4);
TH1F *h_TRG2_totalcharge    =  new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG2_mu_dz0_01      =  new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG2_mu_dz0_12      =  new TH1F(Form("h_TRG%s_mu_dz0_12",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG2_mu_dz0_02      =  new TH1F(Form("h_TRG%s_mu_dz0_02",tagname.c_str() ), "", 80, 0, 20);
TH2F *h_TRG2_respt          =  new TH2F(Form("h_TRG%s_respt"    ,tagname.c_str() ), "", 50, 0, 50, 50, 0, 50);

int trigger2=0;
//int trigger2plus=0;
for (uint i=0; i<GenMuons.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = PureTriggerSelector(Tracks.at(i), mwindow, drwindow, dzwindow);
  if(triggered.size() == 3)
  {
     trigger2++;
     h_TRG2_gen_mu_p_0->Fill(   GenMuons.at(i).at(0).P());
     h_TRG2_gen_mu_pt_0->Fill(  GenMuons.at(i).at(0).Pt());
     h_TRG2_gen_mu_eta_0->Fill( abs(GenMuons.at(i).at(0).Eta())  );
     h_TRG2_gen_mu_p_1->Fill(   GenMuons.at(i).at(1).P());
     h_TRG2_gen_mu_pt_1->Fill(  GenMuons.at(i).at(1).Pt());
     h_TRG2_gen_mu_eta_1->Fill( abs(GenMuons.at(i).at(1).Eta()) );
     h_TRG2_gen_mu_p_2->Fill(   GenMuons.at(i).at(2).P());
     h_TRG2_gen_mu_pt_2->Fill(  GenMuons.at(i).at(2).Pt());
     h_TRG2_gen_mu_eta_2->Fill( abs(GenMuons.at(i).at(2).Eta() ) );
     h_TRG2_gen_mu_p_min->Fill( std::min( { GenMuons.at(i).at(0).P(), GenMuons.at(i).at(1).P(), GenMuons.at(i).at(2).P() }));
     h_TRG2_gen_mu_eta_max->Fill( std::max( { abs(GenMuons.at(i).at(0).Eta()), abs(GenMuons.at(i).at(1).Eta()), abs(GenMuons.at(i).at(2).Eta()) })); 
     h_TRG2_invmass->Fill( (get<0>(triggered.at(0))+get<0>(triggered.at(1))+get<0>(triggered.at(2))).M() );
     h_TRG2_mu_p_0->Fill(   get<0>(triggered.at(0)).P());
     h_TRG2_mu_pt_0->Fill(  get<0>(triggered.at(0)).Pt());
     h_TRG2_mu_eta_0->Fill( abs(get<0>(triggered.at(0)).Eta())  );
     h_TRG2_mu_p_1->Fill(   get<0>(triggered.at(1)).P());
     h_TRG2_mu_pt_1->Fill(  get<0>(triggered.at(1)).Pt());
     h_TRG2_mu_eta_1->Fill( abs(get<0>(triggered.at(1)).Eta())  );
     h_TRG2_mu_p_2->Fill(   get<0>(triggered.at(2)).P() ); 
     h_TRG2_mu_pt_2->Fill(  get<0>(triggered.at(2)).Pt());
     h_TRG2_mu_eta_2->Fill( abs(get<0>(triggered.at(2)).Eta())  );
     h_TRG2_variable->Fill(  get<0>(triggered.at(2)).Pt());
     h_TRG2_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
     h_TRG2_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
     h_TRG2_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
     h_TRG2_mu_dz0_12->Fill( abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) )  );
     h_TRG2_mu_dz0_02->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )  );     
     h_TRG2_respt->Fill( (GenMuons.at(i).at(0)+GenMuons.at(i).at(1)+GenMuons.at(i).at(2)).Pt(), (get<0>(triggered.at(0)) + get<0>(triggered.at(1)) + get<0>(triggered.at(2)) ).Pt() );
  }
  else
  {
     h_TRG2_variable->Fill(-1);  
  }

}

MakeRatePlot(h_TRG2_variable,1,"roc");

//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Efficiency for Trigger %s = ",tagname.c_str() )<<(float)trigger2/GenMuons.size()<<endl;
}


void MixedME0TriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,
  std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12,
  std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float dzwindow, float drwindow, float bendcut, string tagname)
{
//Book Histogram for the reference
//MAKE BINS
Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0};
int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
//Make Histograms for the reference
TH1F *h_TRG3_gen_mu_p_0     =  new TH1F(Form("h_TRG%s_gen_mu_p_0",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG3_gen_mu_pt_0    =  new TH1F(Form("h_TRG%s_gen_mu_pt_0",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG3_gen_mu_eta_0   =  new TH1F(Form("h_TRG%s_gen_mu_eta_0",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG3_gen_mu_p_1     =  new TH1F(Form("h_TRG%s_gen_mu_p_1",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG3_gen_mu_pt_1    =  new TH1F(Form("h_TRG%s_gen_mu_pt_1",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG3_gen_mu_eta_1   =  new TH1F(Form("h_TRG%s_gen_mu_eta_1",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG3_gen_mu_p_2     =  new TH1F(Form("h_TRG%s_gen_mu_p_2",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG3_gen_mu_pt_2    =  new TH1F(Form("h_TRG%s_gen_mu_pt_2",tagname.c_str() ), "",nptbins,ptbins);
TH1F *h_TRG3_gen_mu_eta_2   =  new TH1F(Form("h_TRG%s_gen_mu_eta_2",tagname.c_str() ),"", netabins, etabins);
TH1F *h_TRG3_gen_mu_p_min   =  new TH1F(Form("h_TRG%s_gen_mu_p_min",tagname.c_str() ), "", npbins,pbins);
TH1F *h_TRG3_gen_mu_eta_max =  new TH1F(Form("h_TRG%s_gen_mu_eta_max",tagname.c_str() ), "", netabins,etabins);
TH1F *h_TRG3_mu_p_0         =  new TH1F(Form("h_TRG%s_mu_p_0",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG3_mu_pt_0        =  new TH1F(Form("h_TRG%s_mu_pt_0",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG3_mu_eta_0       =  new TH1F(Form("h_TRG%s_mu_eta_0",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG3_mu_p_1         =  new TH1F(Form("h_TRG%s_mu_p_1",tagname.c_str() ), "", 40, 0, 40);
TH1F *h_TRG3_mu_pt_1        =  new TH1F(Form("h_TRG%s_mu_pt_1",tagname.c_str() ), "",40, 0, 20);
TH1F *h_TRG3_mu_eta_1       =  new TH1F(Form("h_TRG%s_mu_eta_1",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG3_mu_eta_2       =  new TH1F(Form("h_TRG%s_mu_eta_2",tagname.c_str() ),"",40, 0, 4);
TH1F *h_TRG3_maxdR          =  new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "", 100, 0,  4);
TH1F *h_TRG3_totalcharge    =  new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG3_mu_dz0_01      =  new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG3_bend           =  new TH1F(Form("h_TRG%s_bend",  tagname.c_str() ),  "", 50,  0, 0.5);
TH1F *h_TRG3_inversebend    =  new TH1F(Form("h_TRG%s_inversebend",  tagname.c_str() ),  "", 100, 0, 100);
TH1F *h_TRG3_variable       =  new TH1F(Form("h_TRG%s_variable",tagname.c_str() ), "",       400, 0, 100);

int trigger3=0;
//int trigger1plus=0;
for (uint i=0; i<GenMuons.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = MixedME0TriggerSelector(Tracks12.at(i),Tracks3.at(i), dzwindow, drwindow, bendcut); 
  if(triggered.size() == 3){
     trigger3++;
     h_TRG3_gen_mu_p_0->Fill(   GenMuons.at(i).at(0).P());
     h_TRG3_gen_mu_pt_0->Fill(  GenMuons.at(i).at(0).Pt());
     h_TRG3_gen_mu_eta_0->Fill( abs(GenMuons.at(i).at(0).Eta()) );
     h_TRG3_gen_mu_p_1->Fill(   GenMuons.at(i).at(1).P());
     h_TRG3_gen_mu_pt_1->Fill(  GenMuons.at(i).at(1).Pt());
     h_TRG3_gen_mu_eta_1->Fill( abs(GenMuons.at(i).at(1).Eta()) );
     h_TRG3_gen_mu_p_2->Fill(   GenMuons.at(i).at(2).P());
     h_TRG3_gen_mu_pt_2->Fill(  GenMuons.at(i).at(2).Pt());
     h_TRG3_gen_mu_eta_2->Fill( abs(GenMuons.at(i).at(2).Eta()) );
     h_TRG3_gen_mu_p_min->Fill( std::min( { GenMuons.at(i).at(0).P(), GenMuons.at(i).at(1).P(), GenMuons.at(i).at(2).P() }));
     h_TRG3_gen_mu_eta_max->Fill( std::max( { abs(GenMuons.at(i).at(0).Eta()), abs(GenMuons.at(i).at(1).Eta()), abs(GenMuons.at(i).at(2).Eta()) })); 
     h_TRG3_mu_p_0->Fill(   get<0>(triggered.at(0)).P());
     h_TRG3_mu_pt_0->Fill(  get<0>(triggered.at(0)).Pt());
     h_TRG3_mu_eta_0->Fill( abs(get<0>(triggered.at(0)).Eta()) );
     h_TRG3_mu_p_1->Fill(   get<0>(triggered.at(1)).P());
     h_TRG3_mu_pt_1->Fill(  get<0>(triggered.at(1)).Pt());
     h_TRG3_mu_eta_1->Fill( abs(get<0>(triggered.at(1)).Eta()) );
     h_TRG3_mu_eta_2->Fill( abs(get<0>(triggered.at(2)).Eta()) );  
     h_TRG3_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
     h_TRG3_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
     h_TRG3_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
     h_TRG3_bend->Fill(  abs(get<2>(triggered.at(2)) )  );  
     h_TRG3_inversebend->Fill( 1 / abs(get<2>(triggered.at(2)) )  ); 
     h_TRG3_variable->Fill(   1 / abs(get<2>(triggered.at(2)) ) );
  }
  else
  {
      h_TRG3_variable->Fill(-1);    
  }
}
MakeRatePlot(h_TRG3_variable,1,"roc");

//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Efficiency for Trigger %s using %.1f bending =",tagname.c_str(),bendcut)<<(float)trigger3/GenMuons.size()<<endl;

}

void MuonTriggerEfficiencies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> EMTFMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
	std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons ,  
	string outputname)
{
cout << "[INFO] MuonTrigger module running ..."<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_muontriggers.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting Trigger Studies"<<endl;
TriggerEfficiencyStudies(GenMuons,EMTFMuons,L1TTMuons,L1TkMuMuons,L1TkMuStubMuons,ME0StubMuons);

tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  
