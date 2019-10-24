#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TChain.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>   
#include "SelectionHelper.h"
#include "GenProperties.h"
#include <string>
#include "CfgParser.h"
using namespace std;

float CalculateAcceptance(std::vector<std::vector<TLorentzVector>> GenMuons, float muonid, float minpt, float maxeta)
{
//Calculate denominator: Number of events passing our selection
int denominator=GenMuons.size();
//Calculate numberator: Number of events passing pt and eta cuts
int numerator=0;	
for (int i=0; i<denominator; ++i)
{	
if (GenMuons.at(i).at(muonid).Pt() > minpt && GenMuons.at(i).at(muonid).Eta() < maxeta){numerator++;}
}
//Calculate acceptance
float acceptance = (float)numerator/(float)denominator;
return acceptance;
}


void CreateAcceptanceMapHistogram(std::vector<std::vector<TLorentzVector>> GenMuons, float muonid,std::vector<float> pts, std::vector<float> etas,const string name)
{
float effij=0;
TH2F *h_tmp = new TH2F(Form("%s",name.c_str()),Form("%s",name.c_str()),pts.size(),0,pts.size(),etas.size()-1,0,etas.size()-1); 
for (uint i=0; i<pts.size(); ++i)
{	
   for (uint j=1; j<etas.size(); ++j)
   {	
   	effij=(float)(((int)(CalculateAcceptance(GenMuons,muonid,pts.at(i), etas.at(j))*100000))/100000.0);
   	h_tmp->Fill(Form("P_{T}>%.1f",pts.at(i)),Form("%.1f<#eta<%.1f",etas.at(0),etas.at(j)),effij);
   }
}

}

void SelectionStudies(std::vector<std::vector<TLorentzVector>> GenMuons, string outputname)
{
//Event Selection Cuts
vector<float> muonptcuts ={2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};
vector<float> muonetacuts={1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8};
//Create Map Histograms 
CreateAcceptanceMapHistogram(GenMuons,0,muonptcuts,muonetacuts,"First_pT_muon_acceptance_map");
CreateAcceptanceMapHistogram(GenMuons,1,muonptcuts,muonetacuts,"Second_pT_muon_acceptance_map");
CreateAcceptanceMapHistogram(GenMuons,2,muonptcuts,muonetacuts,"Third_pT_muon_acceptance_map");

}

void PropertiesStudies(std::vector<std::vector<TLorentzVector>> gentaus, std::vector<std::vector<TLorentzVector>> genmuons, string outputname)
{
//Book Histograms
//KINEMATIC DISTRIBUTIONS (GEN)
TH1F *h_gen_ev       =  new TH1F("h_gen_ev",        "", 10, 0, 10);
TH1F *h_gen_tau      =  new TH1F("h_gen_tau",       "", 20, 0, 20);
TH1F *h_gen_tau_p_0  =  new TH1F("h_gen_tau_p_0", "", 100, 0,100);
TH1F *h_gen_tau_pt_0 =  new TH1F("h_gen_tau_pt_0", "", 100, 0, 40);
TH1F *h_gen_tau_eta_0=  new TH1F("h_gen_tau_eta_0", "",100, -5, 5);
TH1F *h_gen_tau_phi_0=  new TH1F("h_gen_tau_phi_0", "", 20,-4,  4);
TH1F *h_gen_mu      =   new TH1F("h_gen_mu", "", 20, 0, 20);
TH1F *h_gen_mu_p_0  =   new TH1F("h_gen_mu_p_0", "", 40, 0, 40);
TH1F *h_gen_mu_p_1  =   new TH1F("h_gen_mu_p_1", "", 40, 0, 40);
TH1F *h_gen_mu_p_2  =   new TH1F("h_gen_mu_p_2", "", 40, 0, 40);
TH1F *h_gen_mu_pt_0 =   new TH1F("h_gen_mu_pt_0", "", 40, 0, 20);
TH1F *h_gen_mu_pt_1 =   new TH1F("h_gen_mu_pt_1", "", 40, 0, 20);
TH1F *h_gen_mu_pt_2 =   new TH1F("h_gen_mu_pt_2", "", 40, 0, 20);
TH1F *h_gen_mu_eta_0=   new TH1F("h_gen_mu_eta_0", "", 40, -4, 4);
TH1F *h_gen_mu_eta_1=   new TH1F("h_gen_mu_eta_1", "", 40, -4, 4);
TH1F *h_gen_mu_eta_2=   new TH1F("h_gen_mu_eta_2", "", 40, -4, 4);
TH1F *h_gen_mu_phi_0=   new TH1F("h_gen_mu_phi_0", "", 20, -4, 4);
TH1F *h_gen_mu_phi_1=   new TH1F("h_gen_mu_phi_1", "", 20, -4, 4);
TH1F *h_gen_mu_phi_2=   new TH1F("h_gen_mu_phi_2", "", 20, -4, 4);
TH1F *h_gen_mu_dR01 =   new TH1F("h_gen_mu_dR01", "", 40, 0, 4);
TH1F *h_gen_mu_dR02 =   new TH1F("h_gen_mu_dR02", "", 40, 0, 4);
TH1F *h_gen_mu_dR12 =   new TH1F("h_gen_mu_dR12", "", 40, 0, 4);
TH1F *h_gen_mu_dRmax=   new TH1F("h_gen_mu_dRmax","",40, 0, 4);
TH1F *h_gen_mu_dEta01 = new TH1F("h_gen_mu_dEta01", "", 40, 0, 4);
TH1F *h_gen_mu_dEta02 = new TH1F("h_gen_mu_dEta02", "", 40, 0, 4);
TH1F *h_gen_mu_dEta12 = new TH1F("h_gen_mu_dEta12", "", 40, 0, 4);
TH1F *h_gen_mu_dEtamax= new TH1F("h_gen_mu_dEtamax","",40, 0, 4);
TH1F *h_gen_mu_dPhi01=  new TH1F("h_gen_mu_dPhi01", "", 40, 0, 4);
TH1F *h_gen_mu_dPhi02=  new TH1F("h_gen_mu_dPhi02", "", 40, 0, 4);
TH1F *h_gen_mu_dPhi12=  new TH1F("h_gen_mu_dPhi12", "", 40, 0, 4);
TH1F *h_gen_mu_dPhimax= new TH1F("h_gen_mu_dPhimax","",40, 0, 4);
TH1F *h_gen_invmass   = new TH1F("h_gen_invmass","",100, 0, 3);
TH1F *h_gen_mu_p_min  = new TH1F("h_gen_mu_p_min", "", 100, 0, 20);
TH1F *h_gen_mu_eta_max= new TH1F("h_gen_mu_eta_max", "", 100, -5, 5);
TH2F *h_gen_mu_dEta01_dPhi01 = new TH2F("h_gen_mu_dEta01_dPhi01", "", 20, 0, 2, 50, -5, 5);
TH2F *h_gen_mu_dEta02_dPhi02 = new TH2F("h_gen_mu_dEta02_dPhi02", "", 20, 0, 2, 50, -5, 5);
TH2F *h_gen_mu_dEta12_dPhi12 = new TH2F("h_gen_mu_dEta12_dPhi12", "", 20, 0, 2, 50, -5, 5);
TH2F *h_gen_mu_p_min_eta_max = new TH2F("h_gen_mu_p_min_eta_max", "", 100, 0, 20, 100, -5, 5);
TH2F *h_gen_mu_pt_0_p_0      = new TH2F("h_gen_mu_pt_0_p_0", "",      50,  0, 50,  50, 0, 50);
TH2F *h_gen_mu_pt_1_p_1      = new TH2F("h_gen_mu_pt_1_p_1", "",      50,  0, 50,  50, 0, 50);
TH2F *h_gen_mu_pt_2_p_2      = new TH2F("h_gen_mu_pt_2_p_2", "",      50,  0, 50,  50, 0, 50);
TH2F *h_gen_mu_pt_0_eta_0    = new TH2F("h_gen_mu_pt_0_eta_0", "",    50,  0, 50,  50, -5, 5);
TH2F *h_gen_mu_pt_1_eta_1    = new TH2F("h_gen_mu_pt_1_eta_1", "",    50,  0, 50,  50, -5, 5);
TH2F *h_gen_mu_pt_2_eta_2    = new TH2F("h_gen_mu_pt_2_eta_2", "",    50,  0, 50,  50, -5, 5);

for (uint j=1; j<gentaus.size(); ++j)
{  

//TH1 HISTOGRAMS
h_gen_ev->Fill(0);
h_gen_tau->Fill(       gentaus.at(j).size());
h_gen_tau_p_0->Fill(   gentaus.at(j).at(0).P());
h_gen_tau_pt_0->Fill(  gentaus.at(j).at(0).Pt());
h_gen_tau_eta_0->Fill( gentaus.at(j).at(0).Eta());
h_gen_tau_phi_0->Fill( gentaus.at(j).at(0).Phi());
h_gen_mu->Fill(genmuons.at(j).size());
h_gen_mu_p_0->Fill( genmuons.at(j).at(0).P());
h_gen_mu_p_1->Fill( genmuons.at(j).at(1).P());
h_gen_mu_p_2->Fill( genmuons.at(j).at(2).P());
h_gen_mu_pt_0->Fill( genmuons.at(j).at(0).Pt());
h_gen_mu_pt_1->Fill( genmuons.at(j).at(1).Pt());
h_gen_mu_pt_2->Fill( genmuons.at(j).at(2).Pt());
h_gen_mu_eta_0->Fill( genmuons.at(j).at(0).Eta());
h_gen_mu_eta_1->Fill( genmuons.at(j).at(1).Eta());
h_gen_mu_eta_2->Fill( genmuons.at(j).at(2).Eta());
h_gen_mu_phi_0->Fill( genmuons.at(j).at(0).Phi());
h_gen_mu_phi_1->Fill( genmuons.at(j).at(1).Phi());
h_gen_mu_phi_2->Fill( genmuons.at(j).at(2).Phi());
h_gen_mu_dR01->Fill( genmuons.at(j).at(0).DeltaR(genmuons.at(j).at(1)));
h_gen_mu_dR02->Fill( genmuons.at(j).at(0).DeltaR(genmuons.at(j).at(2)));
h_gen_mu_dR12->Fill( genmuons.at(j).at(1).DeltaR(genmuons.at(j).at(2)));
h_gen_mu_dRmax->Fill(std::max( {genmuons.at(j).at(0).DeltaR(genmuons.at(j).at(1)), genmuons.at(j).at(0).DeltaR(genmuons.at(j).at(2)) , genmuons.at(j).at(1).DeltaR(genmuons.at(j).at(2)  )})   );
h_gen_mu_dPhi01->Fill( genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(1)  ));
h_gen_mu_dPhi02->Fill( genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(2)  ));
h_gen_mu_dPhi12->Fill( genmuons.at(j).at(1).DeltaPhi(genmuons.at(j).at(2)  ));
h_gen_mu_dPhimax->Fill(std::max( {genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(1)), genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(2)), genmuons.at(j).at(1).DeltaPhi(genmuons.at(j).at(2))})   );
h_gen_mu_dEta01->Fill( abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(1).Eta() ) );
h_gen_mu_dEta02->Fill( abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(2).Eta() ) );
h_gen_mu_dEta12->Fill( abs( genmuons.at(j).at(1).Eta() - genmuons.at(j).at(2).Eta() ) );
h_gen_mu_dEtamax->Fill(std::max( {abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(1).Eta() ) , abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(2).Eta() ) , abs( genmuons.at(j).at(1).Eta() - genmuons.at(j).at(2).Eta() )   })   );
h_gen_invmass->Fill( (genmuons.at(j).at(0) + genmuons.at(j).at(1) + genmuons.at(j).at(2) ).M()  );
h_gen_mu_p_min->Fill( std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }));
h_gen_mu_eta_max->Fill( std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }));
h_gen_mu_dEta01_dPhi01->Fill( abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(1).Eta() ), genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(1)) );
h_gen_mu_dEta02_dPhi02->Fill( abs( genmuons.at(j).at(0).Eta() - genmuons.at(j).at(2).Eta() ), genmuons.at(j).at(0).DeltaPhi(genmuons.at(j).at(2)) );
h_gen_mu_dEta12_dPhi12->Fill( abs( genmuons.at(j).at(1).Eta() - genmuons.at(j).at(2).Eta() ), genmuons.at(j).at(1).DeltaPhi(genmuons.at(j).at(2)) );

//TH2 HISTOGRAMS
h_gen_mu_p_min_eta_max->Fill(std::min( { genmuons.at(j).at(0).P(), genmuons.at(j).at(1).P(), genmuons.at(j).at(2).P() }),std::max( { abs(genmuons.at(j).at(0).Eta()), abs(genmuons.at(j).at(1).Eta()), abs(genmuons.at(j).at(2).Eta()) }) );
h_gen_mu_pt_0_p_0->Fill(genmuons.at(j).at(0).Pt(), genmuons.at(j).at(0).P() );
h_gen_mu_pt_1_p_1->Fill(genmuons.at(j).at(1).Pt(), genmuons.at(j).at(1).P() );
h_gen_mu_pt_2_p_2->Fill(genmuons.at(j).at(2).Pt(), genmuons.at(j).at(2).P() );
h_gen_mu_pt_0_eta_0->Fill(genmuons.at(j).at(0).Pt(), genmuons.at(j).at(0).Eta() );
h_gen_mu_pt_1_eta_1->Fill(genmuons.at(j).at(1).Pt(), genmuons.at(j).at(1).Eta() );
h_gen_mu_pt_2_eta_2->Fill(genmuons.at(j).at(2).Pt(), genmuons.at(j).at(2).Eta() );

}

}

void GenProperties(std::vector<std::vector<TLorentzVector>> GenTaus,std::vector<std::vector<TLorentzVector>> GenMuons, string outputname)
{
cout << "[INFO] GenProperties module running ..."<<endl;
cout << " ... Getting only generated taus & muons passing our preselection"<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_genproperties.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting Properties Studies Map      "<<endl;
PropertiesStudies(GenTaus,GenMuons,outputname);
cout << " ... Getting Selection Studies Histograms"<<endl;
SelectionStudies(GenMuons,outputname); 
tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  
