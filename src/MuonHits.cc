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
#include "MuonHits.h"
#include <string>
#include "CfgParser.h"

using namespace std;

void HitPropertiesStudies(std::vector<std::vector<std::array<float,11>>> Hits)
{
//Station S1 histograms
GetHitOccupancies(Hits,1,1,"CSCS1_hits");
GetHitOccupancies(Hits,2,1,"RPCS1_hits");
//GetHitOccupancies(Hits,1,3,"GEMS1_hits");
//GetHitOccupancies(Hits,1,4,"ME0S1_hits");
//Station S2 histograms
GetHitOccupancies(Hits,1,2,"CSCS2_hits");
GetHitOccupancies(Hits,2,2,"RPCS2_hits");
//GetHitOccupancies(Hits,2,3,"GEMS2_hits");
//GetHitOccupancies(Hits,2,4,"ME0S2_hits");
//Station S3 histograms
GetHitOccupancies(Hits,1,3,"CSCS3_hits");
GetHitOccupancies(Hits,2,3,"RPCS3_hits");
//Station S4 histograms
GetHitOccupancies(Hits,1,4,"CSCS4_hits");
GetHitOccupancies(Hits,2,4,"RPCS4_hits");
}


void GetHitOccupancies(std::vector<std::vector<std::array<float,11>>> Hits, int type, int station, string tagname){

//Occupancies by type and station 
TH1F *h_type_station_hits_multiplicity  =  new TH1F(Form("h_%s_multiplicity",tagname.c_str()), "",20,0, 20);
TH1F *h_type_station_hits_phi  =   new TH1F(Form("h_%s_phi",tagname.c_str()), "",  40,-4, 4);
TH1F *h_type_station_hits_eta  =   new TH1F(Form("h_%s_eta",tagname.c_str()), "",  40, 0, 4);
TH1F *h_type_station_hits_r    =   new TH1F(Form("h_%s_r",tagname.c_str()), "", 700,0, 700);
TH1F *h_type_station_hits_z    =   new TH1F(Form("h_%s_z",tagname.c_str()), "", 700,500,1200);
TH1F *h_type_station_hits_x    =   new TH1F(Form("h_%s_x",tagname.c_str()), "", 700,-700, 700);
TH1F *h_type_station_hits_y    =   new TH1F(Form("h_%s_y",tagname.c_str()), "", 700,-700, 700);
TH2F *h_type_station_hits_x_y  =   new TH2F(Form("h_%s_x_y",tagname.c_str()), "",350,-700,700,350,-700,700);
TH2F *h_type_station_hits_z_r  =   new TH2F(Form("h_%s_z_r",tagname.c_str()), "",350,500,1200,350,0,700);
TH2F *h_type_station_hits_z_y  =   new TH2F(Form("h_%s_z_y",tagname.c_str()), "",350,500,1200,350,-700,700);

//Loop over events
for (uint i=1; i<Hits.size(); ++i)
{
   int multiplicity=0; 
   //Loop over hits
   for (uint j=1; j<Hits.at(i).size(); ++j)
   {   
    if(int(Hits.at(i).at(j)[6])==type &&  int(Hits.at(i).at(j)[7])==station){
    multiplicity++;
    h_type_station_hits_x->Fill(  Hits.at(i).at(j)[0] ); 
    h_type_station_hits_y->Fill(  Hits.at(i).at(j)[1] );  
    h_type_station_hits_z->Fill(  Hits.at(i).at(j)[2] ); 
    h_type_station_hits_r->Fill(  Hits.at(i).at(j)[3] );     
    h_type_station_hits_eta->Fill(Hits.at(i).at(j)[4] ); 
    h_type_station_hits_phi->Fill(Hits.at(i).at(j)[5] ); 
    h_type_station_hits_x_y->Fill(Hits.at(i).at(j)[0],Hits.at(i).at(j)[1]);
    h_type_station_hits_z_r->Fill(Hits.at(i).at(j)[2],Hits.at(i).at(j)[3]);
    h_type_station_hits_z_y->Fill(Hits.at(i).at(j)[2],Hits.at(i).at(j)[1]);
    }
   }
   h_type_station_hits_multiplicity->Fill(multiplicity); 
}

}

void MuonHits(std::vector<std::vector<std::array<float,11>>> Hits, string outputname)
{
cout << "[INFO] Muon module running ..."<<endl;
cout << " ... Getting only reconstructed hits passing our preselection"<<endl;
//create output file
TFile *tFile=new TFile( Form("histograms/%s_hits.root",outputname.c_str()), "RECREATE");
//-------------WRITE AND SAVE OUTPUT FILE ---STARTS------------------------------------
//Get generated muons passing event selectiorn
cout << " ... Getting  Hit Properties Studies "<<endl;
HitPropertiesStudies(Hits);
tFile->Write();
tFile->Close();
//-------------WRITE AND SAVE OUTPUT FILE ---ENDS------------------------------------
}  
