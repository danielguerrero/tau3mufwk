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
#include "MuonTrackProperties.h"
#include <string>
#include "CfgParser.h"

using namespace std;
 
void PropertiesStudies(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> tracks, string tagname)
{
//Book Histograms
//KINEMATIC DISTRIBUTIONS (0,1,2,3,more tracks)
TH1F *h_track_mu_n         =   new TH1F(Form("h_%s_mu_n",tagname.c_str()), "",   100, 0, 100);
//KINEMATIC DISTRIBUTIONS (at least 1 track)
TH1F *h_track_mu_pt    = new TH1F(Form("h_%s_mu_pt",tagname.c_str()), "",  40, 0, 20);
TH1F *h_track_mu_eta   = new TH1F(Form("h_%s_mu_eta",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_phi   = new TH1F(Form("h_%s_mu_phi",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_pt_0  = new TH1F(Form("h_%s_mu_pt_0",tagname.c_str()), "",  100, 0, 100);
TH1F *h_track_mu_eta_0 = new TH1F(Form("h_%s_mu_eta_0",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_phi_0 = new TH1F(Form("h_%s_mu_phi_0",tagname.c_str()), "", 20, -4, 4);
//KINEMATIC DISTRIBUTIONS (at least two tracks)
TH1F *h_track_mu_2trk_pt   =   new TH1F(Form("h_%s_mu_2trk_pt",tagname.c_str()), "", 40, 0, 40);
TH1F *h_track_mu_2trk_invmass= new TH1F(Form("h_%s_mu_2trk_invmass",tagname.c_str()), "", 100, 0, 3);
TH1F *h_track_mu_2trk_pt_0 =   new TH1F(Form("h_%s_mu_2trk_pt_0",tagname.c_str()), "", 40, 0, 20);
TH1F *h_track_mu_2trk_pt_1 =   new TH1F(Form("h_%s_mu_2trk_pt_1",tagname.c_str()), "", 40, 0, 20);
TH1F *h_track_mu_2trk_eta_0=   new TH1F(Form("h_%s_mu_2trk_eta_0",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_2trk_eta_1=   new TH1F(Form("h_%s_mu_2trk_eta_1",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_2trk_phi_0=   new TH1F(Form("h_%s_mu_2trk_phi_0",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_2trk_phi_1=   new TH1F(Form("h_%s_mu_2trk_phi_1",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_2trk_dR01 =   new TH1F(Form("h_%s_mu_2trk_dR01",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_2trk_dEta01 = new TH1F(Form("h_%s_mu_2trk_dEta01",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_2trk_dPhi01=  new TH1F(Form("h_%s_mu_2trk_dPhi01",tagname.c_str()), "", 60, 0, 6);
//KINEMATIC DISTRIBUTIONS (at least three tracks)
TH1F *h_track_mu_3trk_pt   =   new TH1F(Form("h_%s_mu_3trk_pt",tagname.c_str()), "", 40, 0, 40);
TH1F *h_track_mu_3trk_invmass= new TH1F(Form("h_%s_mu_3trk_invmass",tagname.c_str()), "", 100, 0, 3);
TH1F *h_track_mu_3trk_pt_0 =   new TH1F(Form("h_%s_mu_3trk_pt_0",tagname.c_str()), "", 40, 0, 20);
TH1F *h_track_mu_3trk_pt_1 =   new TH1F(Form("h_%s_mu_3trk_pt_1",tagname.c_str()), "", 40, 0, 20);
TH1F *h_track_mu_3trk_pt_2 =   new TH1F(Form("h_%s_mu_3trk_pt_2",tagname.c_str()), "", 40, 0, 20);
TH1F *h_track_mu_3trk_eta_0=   new TH1F(Form("h_%s_mu_3trk_eta_0",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_3trk_eta_1=   new TH1F(Form("h_%s_mu_3trk_eta_1",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_3trk_eta_2=   new TH1F(Form("h_%s_mu_3trk_eta_2",tagname.c_str()), "", 30, 0, 3);
TH1F *h_track_mu_3trk_phi_0=   new TH1F(Form("h_%s_mu_3trk_phi_0",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_3trk_phi_1=   new TH1F(Form("h_%s_mu_3trk_phi_1",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_3trk_phi_2=   new TH1F(Form("h_%s_mu_3trk_phi_2",tagname.c_str()), "", 20, -4, 4);
TH1F *h_track_mu_3trk_dR01 =   new TH1F(Form("h_%s_mu_3trk_dR01",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dR02 =   new TH1F(Form("h_%s_mu_3trk_dR02",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dR12 =   new TH1F(Form("h_%s_mu_3trk_dR12",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dRmax=   new TH1F(Form("h_%s_mu_3trk_dRmax",tagname.c_str()),"",60, 0, 6);
TH1F *h_track_mu_3trk_dEta01 = new TH1F(Form("h_%s_mu_3trk_dEta01",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dEta02 = new TH1F(Form("h_%s_mu_3trk_dEta02",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dEta12 = new TH1F(Form("h_%s_mu_3trk_dEta12",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dEtamax= new TH1F(Form("h_%s_mu_3trk_dEtamax",tagname.c_str()),"",60, 0, 6);
TH1F *h_track_mu_3trk_dPhi01=  new TH1F(Form("h_%s_mu_3trk_dPhi01",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dPhi02=  new TH1F(Form("h_%s_mu_3trk_dPhi02",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dPhi12=  new TH1F(Form("h_%s_mu_3trk_dPhi12",tagname.c_str()), "", 60, 0, 6);
TH1F *h_track_mu_3trk_dPhimax= new TH1F(Form("h_%s_mu_3trk_dPhimax",tagname.c_str()),"",60, 0, 6);
TH1F *h_track_mu_3trk_eta_max= new TH1F(Form("h_%s_mu_3trk_eta_max",tagname.c_str()), "", 100, 0, 10);
//ME0 Histograms
TH1F *h_track_mu_bend_0                = new TH1F( Form("h_%s_mu_bend_0                 ",tagname.c_str()),"", 100, -1, 1);
TH1F *h_track_mu_2trk_bend_0           = new TH1F( Form("h_%s_mu_2trk_bend_0            ",tagname.c_str()),"", 100, -1, 1);
TH1F *h_track_mu_2trk_bend_1           = new TH1F( Form("h_%s_mu_2trk_bend_1            ",tagname.c_str()),"", 100, -1, 1);
TH1F *h_track_mu_3trk_bend_0           = new TH1F( Form("h_%s_mu_3trk_bend_0            ",tagname.c_str()),"", 100, -1, 1);
TH1F *h_track_mu_3trk_bend_1           = new TH1F( Form("h_%s_mu_3trk_bend_1            ",tagname.c_str()),"", 100, -1, 1);
TH1F *h_track_mu_3trk_bend_2           = new TH1F( Form("h_%s_mu_3trk_bend_2            ",tagname.c_str()),"", 100, -1, 1);

TH1F *h_track_mu_charge_0              = new TH1F( Form("h_%s_mu_charge_0               ",tagname.c_str()),"", 4, -2 ,2);
TH1F *h_track_mu_2trk_charge_0         = new TH1F( Form("h_%s_mu_2trk_charge_0          ",tagname.c_str()),"", 4, -2 ,2);
TH1F *h_track_mu_2trk_charge_1         = new TH1F( Form("h_%s_mu_2trk_charge_1          ",tagname.c_str()),"", 4, -2 ,2);
TH1F *h_track_mu_3trk_charge_0         = new TH1F( Form("h_%s_mu_3trk_charge_0          ",tagname.c_str()),"", 4, -2 ,2);
TH1F *h_track_mu_3trk_charge_1         = new TH1F( Form("h_%s_mu_3trk_charge_1          ",tagname.c_str()),"", 4, -2 ,2);
TH1F *h_track_mu_3trk_charge_2         = new TH1F( Form("h_%s_mu_3trk_charge_2          ",tagname.c_str()),"", 4, -2 ,2);

TH1F *h_track_mu_chargebend_0          = new TH1F( Form("h_%s_mu_chargebend_0           ",tagname.c_str()),"", 100,  0 ,1);    
TH1F *h_track_mu_2trk_chargebend_0     = new TH1F( Form("h_%s_mu_2trk_chargebend_0      ",tagname.c_str()),"", 100,  0 ,1);
TH1F *h_track_mu_2trk_chargebend_1     = new TH1F( Form("h_%s_mu_2trk_chargebend_1      ",tagname.c_str()),"", 100,  0 ,1);
TH1F *h_track_mu_3trk_chargebend_0     = new TH1F( Form("h_%s_mu_3trk_chargebend_0      ",tagname.c_str()),"", 100,  0 ,1);
TH1F *h_track_mu_3trk_chargebend_1     = new TH1F( Form("h_%s_mu_3trk_chargebend_1      ",tagname.c_str()),"", 100,  0 ,1);
TH1F *h_track_mu_3trk_chargebend_2     = new TH1F( Form("h_%s_mu_3trk_chargebend_2      ",tagname.c_str()),"", 100,  0 ,1);

TH1F *h_track_mu_chargeoverbend_0      = new TH1F( Form("h_%s_mu_chargeoverbend_0       ",tagname.c_str()),"", 100, -1,1);
TH1F *h_track_mu_2trk_chargeoverbend_0 = new TH1F( Form("h_%s_mu_2trk_chargeoverbend_0  ",tagname.c_str()),"", 100, -1,1);
TH1F *h_track_mu_2trk_chargeoverbend_1 = new TH1F( Form("h_%s_mu_2trk_chargeoverbend_1  ",tagname.c_str()),"", 100, -1,1);
TH1F *h_track_mu_3trk_chargeoverbend_0 = new TH1F( Form("h_%s_mu_3trk_chargeoverbend_0  ",tagname.c_str()),"", 100, -1,1);
TH1F *h_track_mu_3trk_chargeoverbend_1 = new TH1F( Form("h_%s_mu_3trk_chargeoverbend_1  ",tagname.c_str()),"", 100, -1,1);
TH1F *h_track_mu_3trk_chargeoverbend_2 = new TH1F( Form("h_%s_mu_3trk_chargeoverbend_2  ",tagname.c_str()),"", 100, -1,1);


for (uint j=0; j<tracks.size(); ++j)
{  
    //Minimal L1T RECONSTRUCTION INFORMATION (No cut on number of tracks)
    h_track_mu_n->Fill( tracks.at(j).size());
    if (tracks.at(j).size()==0) continue;
    //Event with 1 or more  tracks
    for(uint i=0; i<tracks.at(j).size(); i++)
    {h_track_mu_pt->Fill(  get<0>(tracks.at(j).at(i)).Pt());h_track_mu_eta->Fill( get<0>(tracks.at(j).at(i)).Eta());h_track_mu_phi->Fill( get<0>(tracks.at(j).at(i)).Phi());}
    h_track_mu_pt_0->Fill( get<0>(tracks.at(j).at(0)).Pt());
    h_track_mu_eta_0->Fill(get<0>(tracks.at(j).at(0)).Eta());
    h_track_mu_phi_0->Fill(get<0>(tracks.at(j).at(0)).Phi());
     //ME0 information
    if(tagname.substr()=="ME0Stub")
    {
         h_track_mu_charge_0        ->Fill( get<1>(tracks.at(j).at(0)));
         h_track_mu_bend_0          ->Fill( get<2>(tracks.at(j).at(0)));
         h_track_mu_chargeoverbend_0->Fill( get<1>(tracks.at(j).at(0)) / get<2>(tracks.at(j).at(0)) );
         h_track_mu_chargebend_0    ->Fill(-get<1>(tracks.at(j).at(0))*get<2>(tracks.at(j).at(0)) );
    }
    //Event with 2 or more  tracks
    if(tracks.at(j).size()>=2)
    {
    h_track_mu_2trk_invmass->Fill( (get<0>(tracks.at(j).at(0))+get<0>(tracks.at(j).at(1))).M() );
    h_track_mu_2trk_pt->Fill( (get<0>(tracks.at(j).at(0))+get<0>(tracks.at(j).at(1))).Pt() );
    h_track_mu_2trk_pt_0->Fill( get<0>(tracks.at(j).at(0)).Pt());
    h_track_mu_2trk_pt_1->Fill( get<0>(tracks.at(j).at(1)).Pt());
    h_track_mu_2trk_eta_0->Fill( get<0>(tracks.at(j).at(0)).Eta());
    h_track_mu_2trk_eta_1->Fill( get<0>(tracks.at(j).at(1)).Eta());
    h_track_mu_2trk_phi_0->Fill( get<0>(tracks.at(j).at(0)).Phi());
    h_track_mu_2trk_phi_1->Fill( get<0>(tracks.at(j).at(1)).Phi());
    h_track_mu_2trk_dR01->Fill( get<0>(tracks.at(j).at(0)).DeltaR(get<0>(tracks.at(j).at(1))));
    h_track_mu_2trk_dPhi01->Fill( get<0>(tracks.at(j).at(0)).DeltaPhi(get<0>(tracks.at(j).at(1))));
    h_track_mu_2trk_dEta01->Fill( abs( get<0>(tracks.at(j).at(0)).Eta() - get<0>(tracks.at(j).at(1)).Eta() ) );
    //ME0 information
    if(tagname.substr()=="ME0Stub")
    {
         h_track_mu_2trk_charge_0        ->Fill( get<1>(tracks.at(j).at(0)));
         h_track_mu_2trk_charge_1        ->Fill( get<1>(tracks.at(j).at(1)));
         h_track_mu_2trk_bend_0          ->Fill( get<2>(tracks.at(j).at(0)));
         h_track_mu_2trk_bend_1          ->Fill( get<2>(tracks.at(j).at(1)));
         h_track_mu_2trk_chargeoverbend_0->Fill( get<1>(tracks.at(j).at(0)) / get<2>(tracks.at(j).at(0)) );
         h_track_mu_2trk_chargeoverbend_1->Fill( get<1>(tracks.at(j).at(1)) / get<2>(tracks.at(j).at(1)) );
         h_track_mu_2trk_chargebend_0    ->Fill(-get<1>(tracks.at(j).at(0))*get<2>(tracks.at(j).at(0)) );
         h_track_mu_2trk_chargebend_1    ->Fill(-get<1>(tracks.at(j).at(1))*get<2>(tracks.at(j).at(1)) );
    }
    }
    //Events with 3 or more tracks
    if(tracks.at(j).size()>=3)
    {
    h_track_mu_3trk_invmass->Fill( (get<0>(tracks.at(j).at(0))+get<0>(tracks.at(j).at(1))+get<0>(tracks.at(j).at(2))).M() );
    h_track_mu_3trk_pt->Fill( (get<0>(tracks.at(j).at(0))+get<0>(tracks.at(j).at(1))+get<0>(tracks.at(j).at(2))).Pt() );
    h_track_mu_3trk_pt_0->Fill( get<0>(tracks.at(j).at(0)).Pt());
    h_track_mu_3trk_pt_1->Fill( get<0>(tracks.at(j).at(1)).Pt());
    h_track_mu_3trk_pt_2->Fill( get<0>(tracks.at(j).at(2)).Pt());
    h_track_mu_3trk_eta_0->Fill( get<0>(tracks.at(j).at(0)).Eta());
    h_track_mu_3trk_eta_1->Fill( get<0>(tracks.at(j).at(1)).Eta());
    h_track_mu_3trk_eta_2->Fill( get<0>(tracks.at(j).at(2)).Eta());
    h_track_mu_3trk_phi_0->Fill( get<0>(tracks.at(j).at(0)).Phi());
    h_track_mu_3trk_phi_1->Fill( get<0>(tracks.at(j).at(1)).Phi());
    h_track_mu_3trk_phi_2->Fill( get<0>(tracks.at(j).at(2)).Phi());
    h_track_mu_3trk_dR01->Fill( get<0>(tracks.at(j).at(0)).DeltaR(get<0>(tracks.at(j).at(1))));
    h_track_mu_3trk_dR02->Fill( get<0>(tracks.at(j).at(0)).DeltaR(get<0>(tracks.at(j).at(2))));
    h_track_mu_3trk_dR12->Fill( get<0>(tracks.at(j).at(1)).DeltaR(get<0>(tracks.at(j).at(2))));
    h_track_mu_3trk_dRmax->Fill(std::max( {get<0>(tracks.at(j).at(0)).DeltaR(get<0>(tracks.at(j).at(1))), get<0>(tracks.at(j).at(0)).DeltaR(get<0>(tracks.at(j).at(2))), get<0>(tracks.at(j).at(1)).DeltaR(get<0>(tracks.at(j).at(2)))})   );
    h_track_mu_3trk_dPhi01->Fill( get<0>(tracks.at(j).at(0)).DeltaPhi(get<0>(tracks.at(j).at(1))));
    h_track_mu_3trk_dPhi02->Fill( get<0>(tracks.at(j).at(0)).DeltaPhi(get<0>(tracks.at(j).at(2))));
    h_track_mu_3trk_dPhi12->Fill( get<0>(tracks.at(j).at(1)).DeltaPhi(get<0>(tracks.at(j).at(2))));
    h_track_mu_3trk_dPhimax->Fill(std::max( {get<0>(tracks.at(j).at(0)).DeltaPhi(get<0>(tracks.at(j).at(1))), get<0>(tracks.at(j).at(0)).DeltaPhi(get<0>(tracks.at(j).at(2))), get<0>(tracks.at(j).at(1)).DeltaPhi(get<0>(tracks.at(j).at(2)))})   );
    h_track_mu_3trk_dEta01->Fill( abs( get<0>(tracks.at(j).at(0)).Eta() - get<0>(tracks.at(j).at(1)).Eta() ) );
    h_track_mu_3trk_dEta02->Fill( abs( get<0>(tracks.at(j).at(0)).Eta() - get<0>(tracks.at(j).at(2)).Eta() ) );
    h_track_mu_3trk_dEta12->Fill( abs( get<0>(tracks.at(j).at(1)).Eta() - get<0>(tracks.at(j).at(2)).Eta() ) );
    h_track_mu_3trk_dEtamax->Fill(std::max( {abs( get<0>(tracks.at(j).at(0)).Eta() - get<0>(tracks.at(j).at(1)).Eta() ) , abs( get<0>(tracks.at(j).at(0)).Eta() - get<0>(tracks.at(j).at(2)).Eta() ) , abs( get<0>(tracks.at(j).at(1)).Eta() - get<0>(tracks.at(j).at(2)).Eta() )   })   );
    h_track_mu_3trk_eta_max->Fill( std::max( { abs(get<0>(tracks.at(j).at(0)).Eta()), abs(get<0>(tracks.at(j).at(1)).Eta()), abs(get<0>(tracks.at(j).at(2)).Eta()) }));
    //ME0 information
    if(tagname.substr()=="ME0Stub")
    {
         h_track_mu_3trk_charge_0->Fill(         get<1>(tracks.at(j).at(0)));
         h_track_mu_3trk_charge_1->Fill(         get<1>(tracks.at(j).at(1)));
         h_track_mu_3trk_charge_2->Fill(         get<1>(tracks.at(j).at(2)));
         h_track_mu_3trk_bend_0->Fill(           get<2>(tracks.at(j).at(0)));
         h_track_mu_3trk_bend_1->Fill(           get<2>(tracks.at(j).at(1)));
         h_track_mu_3trk_bend_2->Fill(           get<2>(tracks.at(j).at(2)));
         h_track_mu_3trk_chargeoverbend_0->Fill( get<1>(tracks.at(j).at(0)) / get<2>(tracks.at(j).at(0)) );
         h_track_mu_3trk_chargeoverbend_1->Fill( get<1>(tracks.at(j).at(1)) / get<2>(tracks.at(j).at(1)) );
         h_track_mu_3trk_chargeoverbend_2->Fill( get<1>(tracks.at(j).at(2)) / get<2>(tracks.at(j).at(2)) );
         h_track_mu_3trk_chargebend_0->Fill(    -get<1>(tracks.at(j).at(0))*get<2>(tracks.at(j).at(0)) );
         h_track_mu_3trk_chargebend_1->Fill(    -get<1>(tracks.at(j).at(1))*get<2>(tracks.at(j).at(1)) );
         h_track_mu_3trk_chargebend_2->Fill(    -get<1>(tracks.at(j).at(2))*get<2>(tracks.at(j).at(2)) );
    }
    }

}

}

void MuonTrackProperties(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons,
    std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons,  
	string outputname)
{
cout << "[INFO] MuonTrackProperties module running ..."<<endl;
cout << " ... Getting only tracks passing our preselection in config file"<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_trackproperties.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting Properties Studies: EMTF"<<endl;
PropertiesStudies(EMTFMuons,"EMTF");
cout << " ... Getting Properties Studies: L1TT"<<endl;
PropertiesStudies(L1TTMuons,"L1TT");
cout << " ... Getting Properties Studies: L1TkMu"<<endl;
PropertiesStudies(L1TkMuMuons,"L1TkMu");
cout << " ... Getting Properties Studies: L1TkMuStub"<<endl;
PropertiesStudies(L1TkMuStubMuons,"L1TkMuStub");
cout << " ... Getting Properties Studies: ME0Stub"<<endl;
PropertiesStudies(ME0StubMuons,"ME0Stub");

tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  
