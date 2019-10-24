#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>   
#include "SelectionHelper.h"
#include "MuonEfficiencies.h"
#include <string>
#include "CfgParser.h"
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>

using namespace std;

void EfficienciesStudies(std::vector<std::vector<TLorentzVector>> genmuons, std::vector<std::vector<TLorentzVector>> tracks, string tagname)
{
//Book Histograms 
//MAKE BINS
Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0};
int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
//Make Histograms
TH1F *h_REF_gen_mu_p_0          =  new TH1F(Form("h_REF_%s_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_0         =  new TH1F(Form("h_REF_%s_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_0        =  new TH1F(Form("h_REF_%s_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
TH1F *h_REF_gen_mu_p_1          =  new TH1F(Form("h_REF_%s_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_1         =  new TH1F(Form("h_REF_%s_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_1        =  new TH1F(Form("h_REF_%s_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
TH1F *h_REF_gen_mu_p_2          =  new TH1F(Form("h_REF_%s_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
TH1F *h_REF_gen_mu_pt_2         =  new TH1F(Form("h_REF_%s_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_REF_gen_mu_eta_2        =  new TH1F(Form("h_REF_%s_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins);
TH1F *h_REF_gen_mu_p_min        =  new TH1F(Form("h_REF_%s_gen_mu_p_min",tagname.c_str()), "", npbins,pbins);
TH1F *h_REF_gen_mu_eta_max      =  new TH1F(Form("h_REF_%s_gen_mu_eta_max",tagname.c_str()), "", netabins,etabins);

TH1F *h_track1_gen_mu_p_0          =  new TH1F(Form("h_%s1_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
TH1F *h_track1_gen_mu_pt_0         =  new TH1F(Form("h_%s1_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track1_gen_mu_eta_0        =  new TH1F(Form("h_%s1_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
TH1F *h_track1_gen_mu_p_1          =  new TH1F(Form("h_%s1_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
TH1F *h_track1_gen_mu_pt_1         =  new TH1F(Form("h_%s1_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track1_gen_mu_eta_1        =  new TH1F(Form("h_%s1_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
TH1F *h_track1_gen_mu_p_2          =  new TH1F(Form("h_%s1_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
TH1F *h_track1_gen_mu_pt_2         =  new TH1F(Form("h_%s1_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track1_gen_mu_eta_2        =  new TH1F(Form("h_%s1_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins);
TH1F *h_track1_gen_mu_p_min        =  new TH1F(Form("h_%s1_gen_mu_p_min",tagname.c_str()), "", npbins,pbins);
TH1F *h_track1_gen_mu_eta_max      =  new TH1F(Form("h_%s1_gen_mu_eta_max",tagname.c_str()), "", netabins,etabins);

TH1F *h_track2_gen_mu_p_0          =  new TH1F(Form("h_%s2_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
TH1F *h_track2_gen_mu_pt_0         =  new TH1F(Form("h_%s2_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track2_gen_mu_eta_0        =  new TH1F(Form("h_%s2_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
TH1F *h_track2_gen_mu_p_1          =  new TH1F(Form("h_%s2_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
TH1F *h_track2_gen_mu_pt_1         =  new TH1F(Form("h_%s2_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track2_gen_mu_eta_1        =  new TH1F(Form("h_%s2_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
TH1F *h_track2_gen_mu_p_2          =  new TH1F(Form("h_%s2_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
TH1F *h_track2_gen_mu_pt_2         =  new TH1F(Form("h_%s2_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track2_gen_mu_eta_2        =  new TH1F(Form("h_%s2_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins);
TH1F *h_track2_gen_mu_p_min        =  new TH1F(Form("h_%s2_gen_mu_p_min",tagname.c_str()), "", npbins,pbins);
TH1F *h_track2_gen_mu_eta_max      =  new TH1F(Form("h_%s2_gen_mu_eta_max",tagname.c_str()), "", netabins,etabins);

TH1F *h_track3_gen_mu_p_0          =  new TH1F(Form("h_%s3_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
TH1F *h_track3_gen_mu_pt_0         =  new TH1F(Form("h_%s3_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track3_gen_mu_eta_0        =  new TH1F(Form("h_%s3_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
TH1F *h_track3_gen_mu_p_1          =  new TH1F(Form("h_%s3_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
TH1F *h_track3_gen_mu_pt_1         =  new TH1F(Form("h_%s3_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track3_gen_mu_eta_1        =  new TH1F(Form("h_%s3_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
TH1F *h_track3_gen_mu_p_2          =  new TH1F(Form("h_%s3_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
TH1F *h_track3_gen_mu_pt_2         =  new TH1F(Form("h_%s3_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
TH1F *h_track3_gen_mu_eta_2        =  new TH1F(Form("h_%s3_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins);
TH1F *h_track3_gen_mu_p_min        =  new TH1F(Form("h_%s3_gen_mu_p_min",tagname.c_str()), "", npbins,pbins);
TH1F *h_track3_gen_mu_eta_max      =  new TH1F(Form("h_%s3_gen_mu_eta_max",tagname.c_str()), "", netabins,etabins);

for (uint j=0; j<genmuons.size(); ++j)
{  
//-------------FILL HISTOGRAMS----STARTS---------------------------------------
h_REF_gen_mu_p_0->Fill(   genmuons.at(j).at(0).P());
h_REF_gen_mu_pt_0->Fill(  genmuons.at(j).at(0).Pt());
h_REF_gen_mu_eta_0->Fill( genmuons.at(j).at(0).Eta());
h_REF_gen_mu_p_1->Fill(   genmuons.at(j).at(1).P());
h_REF_gen_mu_pt_1->Fill(  genmuons.at(j).at(1).Pt());
h_REF_gen_mu_eta_1->Fill( genmuons.at(j).at(1).Eta());
h_REF_gen_mu_p_2->Fill(   genmuons.at(j).at(2).P());
h_REF_gen_mu_pt_2->Fill(  genmuons.at(j).at(2).Pt());
h_REF_gen_mu_eta_2->Fill( genmuons.at(j).at(2).Eta());
h_REF_gen_mu_p_min->Fill( std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }));
h_REF_gen_mu_eta_max->Fill( std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }));


if(tracks.at(j).size()>=1) 
{
h_track1_gen_mu_p_0->Fill(   genmuons.at(j).at(0).P());
h_track1_gen_mu_pt_0->Fill(  genmuons.at(j).at(0).Pt());
h_track1_gen_mu_eta_0->Fill( genmuons.at(j).at(0).Eta());
h_track1_gen_mu_p_1->Fill(   genmuons.at(j).at(1).P());
h_track1_gen_mu_pt_1->Fill(  genmuons.at(j).at(1).Pt());
h_track1_gen_mu_eta_1->Fill( genmuons.at(j).at(1).Eta());
h_track1_gen_mu_p_2->Fill(   genmuons.at(j).at(2).P());
h_track1_gen_mu_pt_2->Fill(  genmuons.at(j).at(2).Pt());
h_track1_gen_mu_eta_2->Fill( genmuons.at(j).at(2).Eta());
h_track1_gen_mu_p_min->Fill( std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }));
h_track1_gen_mu_eta_max->Fill( std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }));
}
if(tracks.at(j).size()>=2) 
{
h_track2_gen_mu_p_0->Fill(   genmuons.at(j).at(0).P());
h_track2_gen_mu_pt_0->Fill(  genmuons.at(j).at(0).Pt());
h_track2_gen_mu_eta_0->Fill( genmuons.at(j).at(0).Eta());
h_track2_gen_mu_p_1->Fill(   genmuons.at(j).at(1).P());
h_track2_gen_mu_pt_1->Fill(  genmuons.at(j).at(1).Pt());
h_track2_gen_mu_eta_1->Fill( genmuons.at(j).at(1).Eta());
h_track2_gen_mu_p_2->Fill(   genmuons.at(j).at(2).P());
h_track2_gen_mu_pt_2->Fill(  genmuons.at(j).at(2).Pt());
h_track2_gen_mu_eta_2->Fill( genmuons.at(j).at(2).Eta());
h_track2_gen_mu_p_min->Fill( std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }));
h_track2_gen_mu_eta_max->Fill( std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }));
}
if(tracks.at(j).size()>=3) 
{
h_track3_gen_mu_p_0->Fill(   genmuons.at(j).at(0).P());
h_track3_gen_mu_pt_0->Fill(  genmuons.at(j).at(0).Pt());
h_track3_gen_mu_eta_0->Fill( genmuons.at(j).at(0).Eta());
h_track3_gen_mu_p_1->Fill(   genmuons.at(j).at(1).P());
h_track3_gen_mu_pt_1->Fill(  genmuons.at(j).at(1).Pt());
h_track3_gen_mu_eta_1->Fill( genmuons.at(j).at(1).Eta());
h_track3_gen_mu_p_2->Fill(   genmuons.at(j).at(2).P());
h_track3_gen_mu_pt_2->Fill(  genmuons.at(j).at(2).Pt());
h_track3_gen_mu_eta_2->Fill( genmuons.at(j).at(2).Eta());
h_track3_gen_mu_p_min->Fill( std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }));
h_track3_gen_mu_eta_max->Fill( std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }));
}
//-------------FILL HISTOGRAMS----ENDS---------------------------------------

}

}


void MuonEfficiencies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<TLorentzVector>> EMTFMuons,
	std::vector<std::vector<TLorentzVector>> L1TTMuons,
	std::vector<std::vector<TLorentzVector>> L1TkMuMuons,
	std::vector<std::vector<TLorentzVector>> L1TkMuStubMuons,string outputname)
{
cout << "[INFO] Muon Trigger Efficiencies module running ..."<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_l1tefficiencies.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting Efficiencies Studies     "<<endl;
cout << " ... Getting: EMTF"<<endl;
EfficienciesStudies(GenMuons,EMTFMuons,"EMTF");
cout << " ... Getting: L1TT"<<endl;
EfficienciesStudies(GenMuons,L1TTMuons,"L1TT");
cout << " ... Getting: L1TkMu"<<endl;
EfficienciesStudies(GenMuons,L1TkMuMuons,"L1TKMU");
cout << " ... Getting: L1TkMuStub"<<endl;
EfficienciesStudies(GenMuons,L1TkMuStubMuons,"L1TKMUSTUB");

tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  