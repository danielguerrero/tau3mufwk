#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
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
#include "CfgParser.h"
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>

using namespace std;

std::vector< std::tuple<std::vector<TLorentzVector>,
std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,
std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,
std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,
std::vector<std::array<float,11>> >> SelectParticlesTracksHits(TTree *tree, CfgParser config, const bool is_signal)
{
//Book variables and set branches address
uint         n_gen_mu;        if(is_signal) tree->SetBranchAddress("n_gen_mu",    &(n_gen_mu));
vector<float>* gen_mu_eta=0;  if(is_signal) tree->SetBranchAddress("gen_mu_eta",  &(gen_mu_eta)); 
vector<float>* gen_mu_pt=0;   if(is_signal) tree->SetBranchAddress("gen_mu_pt",   &(gen_mu_pt));
vector<float>* gen_mu_phi=0;  if(is_signal) tree->SetBranchAddress("gen_mu_phi",  &(gen_mu_phi));
vector<float>* gen_mu_e=0;    if(is_signal) tree->SetBranchAddress("gen_mu_e",    &(gen_mu_e)); 
uint         n_gen_tau;       if(is_signal) tree->SetBranchAddress("n_gen_tau",    &(n_gen_tau));
vector<float>* gen_tau_eta=0; if(is_signal) tree->SetBranchAddress("gen_tau_eta",  &(gen_tau_eta));  
vector<float>* gen_tau_pt=0;  if(is_signal) tree->SetBranchAddress("gen_tau_pt",   &(gen_tau_pt));
vector<float>* gen_tau_phi=0; if(is_signal) tree->SetBranchAddress("gen_tau_phi",  &(gen_tau_phi));
vector<float>* gen_tau_e=0;   if(is_signal) tree->SetBranchAddress("gen_tau_e",    &(gen_tau_e));
uint         n_EMTF_mu;           tree->SetBranchAddress("n_EMTF_mu",    &(n_EMTF_mu));  
vector<float>* EMTF_mu_eta=0;     tree->SetBranchAddress("EMTF_mu_eta",  &(EMTF_mu_eta));  
vector<float>* EMTF_mu_pt=0;      tree->SetBranchAddress("EMTF_mu_pt",   &(EMTF_mu_pt));
vector<float>* EMTF_mu_pt_xml=0;  tree->SetBranchAddress("EMTF_mu_pt_xml",   &(EMTF_mu_pt_xml));
vector<float>* EMTF_mu_phi=0;     tree->SetBranchAddress("EMTF_mu_phi",  &(EMTF_mu_phi));
vector<int>* EMTF_mu_charge=0;    tree->SetBranchAddress("EMTF_mu_charge",  &(EMTF_mu_charge));
vector<float>* EMTF_mu_chi2=0;    vector<float>* EMTF_mu_z=0; vector<int>* EMTF_mu_nstubs =0;
uint         n_L1TT_trk;          tree->SetBranchAddress("n_L1TT_trk",    &(n_L1TT_trk));  
vector<float>* L1TT_trk_eta=0;    tree->SetBranchAddress("L1TT_trk_eta",  &(L1TT_trk_eta));  
vector<float>* L1TT_trk_pt=0;     tree->SetBranchAddress("L1TT_trk_pt",   &(L1TT_trk_pt));
vector<float>* L1TT_trk_phi=0;    tree->SetBranchAddress("L1TT_trk_phi",  &(L1TT_trk_phi));
vector<float>* L1TT_trk_chi2=0;   tree->SetBranchAddress("L1TT_trk_chi2",  &(L1TT_trk_chi2));
vector<int>*   L1TT_trk_nstubs=0; tree->SetBranchAddress("L1TT_trk_nstubs",  &(L1TT_trk_nstubs));
vector<int>*   L1TT_trk_charge=0; tree->SetBranchAddress("L1TT_trk_charge",  &(L1TT_trk_charge));
vector<float>* L1TT_trk_z=0;      tree->SetBranchAddress("L1TT_trk_z",     &(L1TT_trk_z));
uint         n_L1_TkMu;           tree->SetBranchAddress("n_L1_TkMu",      &(n_L1_TkMu));  
vector<float>* L1_TkMu_eta=0;     tree->SetBranchAddress("L1_TkMu_eta",    &(L1_TkMu_eta));  
vector<float>* L1_TkMu_pt=0;      tree->SetBranchAddress("L1_TkMu_pt",     &(L1_TkMu_pt));
vector<float>* L1_TkMu_phi=0;     tree->SetBranchAddress("L1_TkMu_phi",    &(L1_TkMu_phi));
vector<int>*   L1_TkMu_charge=0;  tree->SetBranchAddress("L1_TkMu_charge", &(L1_TkMu_charge));
vector<float>* L1_TkMu_z=0;       tree->SetBranchAddress("L1_TkMu_z",  &(L1_TkMu_z));
vector<float>* L1_TkMu_chi2=0;vector<int>* L1_TkMu_nstubs=0;
uint         n_L1_TkMuStub;          tree->SetBranchAddress("n_L1_TkMuStub",      &(n_L1_TkMuStub));  
vector<float>* L1_TkMuStub_eta=0;    tree->SetBranchAddress("L1_TkMuStub_eta",    &(L1_TkMuStub_eta));  
vector<float>* L1_TkMuStub_pt=0;     tree->SetBranchAddress("L1_TkMuStub_pt",     &(L1_TkMuStub_pt));
vector<float>* L1_TkMuStub_phi=0;    tree->SetBranchAddress("L1_TkMuStub_phi",    &(L1_TkMuStub_phi));
vector<int>*   L1_TkMuStub_charge=0; tree->SetBranchAddress("L1_TkMuStub_charge", &(L1_TkMuStub_charge));
vector<float>* L1_TkMuStub_z=0;      tree->SetBranchAddress("L1_TkMuStub_z",      &(L1_TkMuStub_z));
vector<float>* L1_TkMuStub_chi2=0; vector<int>* L1_TkMuStub_nstubs=0;
uint         n_mu_hit;     tree->SetBranchAddress("n_mu_hit",    &(n_mu_hit));
vector<float>* mu_hit_sim_eta=0;    tree->SetBranchAddress("mu_hit_sim_eta",  &(mu_hit_sim_eta));  
vector<float>* mu_hit_sim_phi=0;    tree->SetBranchAddress("mu_hit_sim_phi",  &(mu_hit_sim_phi));
vector<float>* mu_hit_sim_r=0;      tree->SetBranchAddress("mu_hit_sim_r",  &(mu_hit_sim_r));  
vector<float>* mu_hit_sim_z=0;      tree->SetBranchAddress("mu_hit_sim_z",  &(mu_hit_sim_z));
vector<int>* mu_hit_type=0;         tree->SetBranchAddress("mu_hit_type",     &(mu_hit_type));
vector<int>* mu_hit_ring=0;         tree->SetBranchAddress("mu_hit_ring",     &(mu_hit_ring));
vector<int>* mu_hit_station=0;      tree->SetBranchAddress("mu_hit_station", &(mu_hit_station));
vector<int>* mu_hit_sector=0;       tree->SetBranchAddress("mu_hit_sector", &(mu_hit_sector));
vector<int>* mu_hit_neighbor=0;     tree->SetBranchAddress("mu_hit_neighbor",     &(mu_hit_neighbor)); 
vector<float>* mu_hit_bend=0;       tree->SetBranchAddress("mu_hit_bend",     &(mu_hit_bend));
vector<float>* mu_hit_sim_theta=0;  tree->SetBranchAddress("mu_hit_sim_theta",     &(mu_hit_sim_theta));
vector<int>* mu_hit_quality=0;      tree->SetBranchAddress("mu_hit_quality",     &(mu_hit_quality));

//-------LOOP OVER EVENTS---STARTS-------------------------------------
cout << "[INFO] Event processing starts now!"<<endl;
std::vector< std::tuple<std::vector<TLorentzVector>,
std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,
std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,
std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::array<float,11>> >> ParticlesTracksHits;

for (int i=0; i<tree->GetEntries(); ++i)
{ 
//GET EVENT INFORMATION   
tree->GetEvent(i); 
if (i % 1000000 == 0) cout << "    ... processing event " << i << endl;
std::vector<std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>>> genparticles;
std::vector<TLorentzVector> gentaus,genmuons;
int endcap=0;
//If sample is Tau3Mu signal add geninformation, otherwise will be empty vectors
if (is_signal)
{
  gentaus =BuildVectorOfTLorentzVector(gen_tau_pt,gen_tau_eta,gen_tau_phi,gen_tau_e); 
  genmuons=BuildVectorOfTLorentzVector(gen_mu_pt,gen_mu_eta,gen_mu_phi,gen_mu_e);
  //-------------GEN-LEVEL--PRESELECTION---STARTS-------------------------------------
  genparticles = SelectParticles(gentaus, genmuons, config);
  if(genparticles.size()==0) continue; 
  if(get<1>(genparticles[0]).at(0).Eta()>0){endcap =1;}
  else{endcap =-1;}  
  //-------------GEN-LEVEL--PRESELECTION----ENDS---------------------------------------
}

//-------------TRACK--SELECTION---STARTS-------------------------------------
std::vector<std::tuple<TLorentzVector,int,float>> emtfs,l1tts,l1tkmus,l1tkmustubs,em0stubs; 
std::vector<std::array<float,11>> hits;   
emtfs                                   = SelectTracks(EMTF_mu_pt_xml,EMTF_mu_eta,EMTF_mu_phi,EMTF_mu_chi2,EMTF_mu_nstubs,EMTF_mu_charge,EMTF_mu_z,true,false,true,endcap,config); 
l1tts                                   = SelectTracks(L1TT_trk_pt,L1TT_trk_eta,L1TT_trk_phi,L1TT_trk_chi2,L1TT_trk_nstubs,L1TT_trk_charge,L1TT_trk_z,false,true,false,endcap, config);
l1tkmus                                 = SelectTracks(L1_TkMu_pt,L1_TkMu_eta,L1_TkMu_phi,L1_TkMu_chi2,L1_TkMu_nstubs,L1_TkMu_charge,L1_TkMu_z,false,false,false,endcap,config);
l1tkmustubs                             = SelectTracks(L1_TkMuStub_pt,L1_TkMuStub_eta,L1_TkMuStub_phi,L1_TkMuStub_chi2,L1_TkMuStub_nstubs,L1_TkMuStub_charge,L1_TkMuStub_z,false,false,false,endcap,config);
em0stubs                                = SelectStubs( mu_hit_sim_eta,mu_hit_sim_phi,mu_hit_sim_theta,mu_hit_bend,mu_hit_quality,mu_hit_type,endcap);
hits                                    = SelectHits(  mu_hit_sim_eta,mu_hit_sim_phi,mu_hit_sim_r,mu_hit_sim_z,mu_hit_type,mu_hit_station,mu_hit_ring,mu_hit_sector,mu_hit_neighbor,endcap,config);
//-------------TRACK--SELECTION----ENDS---------------------------------------

//Fill Vector of Tuples
if(is_signal)
  ParticlesTracksHits.push_back(make_tuple(get<0>(genparticles[0]), get<1>(genparticles[0]), emtfs, l1tts, l1tkmus, l1tkmustubs, em0stubs, hits));
else
  ParticlesTracksHits.push_back(make_tuple(gentaus, genmuons, emtfs, l1tts, l1tkmus, l1tkmustubs, em0stubs, hits));
}
//-------LOOP OVER EVENTS---STARTS-------------------------------------


//-------RETURN GEN-MUONS----------------------------------------------s
return ParticlesTracksHits;
} 

std::vector<std::array<float,11>> SelectHits(vector<float>* eta, vector<float>* phi, vector<float>* r_xy, vector<float>* z, vector<int>* type, vector<int>* station, vector<int>* ring, vector<int>* sector, vector<int>* neighbor, int endcap, CfgParser config)
{
std::vector<std::array<float,11>> output;
float HITMINETA=stof(config.readStringOpt("Hits::HITMINETA"));
float HITMAXETA=stof(config.readStringOpt("Hits::HITMAXETA"));
for(uint i=0; i<eta->size(); i++)
{
  if(endcap>0)
  {
      if(eta->at(i)>HITMAXETA || eta->at(i)<HITMINETA ) continue; 
  }
  else if(endcap<0)
  {
      if(eta->at(i)<-HITMAXETA || eta->at(i)>-HITMINETA ) continue;    
  }
  else
  {
  //do nothing
  }

  if(neighbor->at(i)==1) continue;
  float x = r_xy->at(i) * TMath::Cos(phi->at(i)* TMath::Pi() / 180);
  float y = r_xy->at(i) * TMath::Sin(phi->at(i)* TMath::Pi() / 180);
  float phiconverted = phi->at(i) * TMath::Pi() / 180;
  std::array<float,11> array; 
  array= {x,y,(float)z->at(i),(float)r_xy->at(i),(float)eta->at(i),phiconverted,(float)type->at(i),(float)station->at(i),(float)ring->at(i),(float)sector->at(i),(float)neighbor->at(i)};
  output.push_back(array);
}

return output;
}

std::vector<std::tuple<float,float,float,int,float>> RemoveOvelapppingTracks(std::vector<std::tuple<float,float,float,int,float>> info)
{
  auto it = info.begin();
  while(it != info.end())
  {
      int found=0;
      TVector3 first,next;
      first.SetPtEtaPhi(get<0>(*it),get<1>(*it),get<2>(*it));
      for(uint k=0;k<info.size();k++)
      {
          next.SetPtEtaPhi(get<0>(info[k]),get<1>(info[k]),get<2>(info[k]));
          if( abs(first.Pt()-next.Pt())<0.1 || abs(first.Eta()-next.Eta())<0.01 || abs(first.DeltaPhi(next))<0.01 )
            {
              found++;
              //if(found>1){cout<<"Pts ("<<first.Pt()<<","<<next.Pt()<<"); Etas ("<<first.Eta()<<","<<next.Eta()<<"); Phis ("<<first.Phi()<<","<<next.Phi()<<")"<<endl;}
            }
      } 
      if(found>1)
      {
      it = info.erase(it);
      continue;    
      }  
      it++;
  }
  return info;
}

std::vector<std::tuple<TLorentzVector,int,float>> SelectStubs(vector<float>* eta, vector<float>* phi, vector<float>* theta, vector<float>* bend, vector<int>* quality, vector<int>* type, int endcap)
{
//Only focused on ME0 (2.0<|eta|<2.8)
float minpt  = 0;
float mineta = 2.0; 
float maxeta = 2.8;
// Note this function is designed to select stubs properties P4, CHARGE, bending(rad)
std::vector<std::tuple<float,float,float,int,float>> tmpoutput,tmpoutput2;
std::vector<std::tuple<TLorentzVector,int,float>> output;
//Make a tuple and remove overlapping track information 
std::vector<std::tuple<float,float,float,int,float>> info;
for(uint k=0;k<type->size();k++)
{
  if( type->at(k) != 4 ) continue; //Only ME0 stubs
  if( quality->at(k) < 6) continue; //Only ME0 Stubs with 5,6 layers
  float pt = ( -0.0039 + (0.134*theta->at(k)) ) / abs( bend->at(k)*TMath::Pi()/180 );
  float charge = -1* TMath::Sign(1, bend->at(k) );
  info.push_back(make_tuple(pt, eta->at(k), (phi->at(k)*TMath::Pi()/180), charge, ( 0.20833333333*(TMath::Pi()/180)*bend->at(k) ) ) );
}
//If no tracks return empty
if(info.size()==0) return output;
//Fill inclusive endcap EMTF muons (for general studies)
if(endcap==0) 
{
    //Fill with positive endcap EMTF muons and ordered them by pt
    uint i=0; 
    while (i<info.size()) 
    {
      if(get<1>(info[i])>mineta && get<1>(info[i])<maxeta && get<0>(info[i])>minpt)
      {            
            tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
      }
      i++;
    }
    std::sort(tmpoutput.rbegin(), tmpoutput.rend());
    //Fill with negative endcap EMTF muons and ordered them by pt
    uint h=0; 
    while (h<info.size()) 
    {
      if(get<1>(info[h])<-mineta && get<1>(info[h])>-maxeta && get<0>(info[h])>minpt)
      {            
            tmpoutput2.push_back(make_tuple(get<0>(info[h]),get<1>(info[h]),get<2>(info[h]),get<3>(info[h]),get<4>(info[h])  ));
      }
      h++;
    }
    std::sort(tmpoutput2.rbegin(), tmpoutput2.rend());
    //Mergue both endcaps
    tmpoutput.insert(tmpoutput.end(),tmpoutput2.begin(),tmpoutput2.end());
}
else
{
    uint i=0; 
    while (i<info.size()) 
    {
      if(endcap>0) //Positive endcap (for signal only)
      {
         if(get<1>(info[i])>mineta && get<1>(info[i])<maxeta && get<0>(info[i])>minpt)
         {            
               tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
         }
         i++;
      }
      else        //Negative endcap (for signal only)
      {
         if(get<1>(info[i])<-mineta && get<1>(info[i])>-maxeta && get<0>(info[i])>minpt)
         {            
               tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
         }
         i++;
      }
    }
    //Order by pt
    std::sort(tmpoutput.rbegin(), tmpoutput.rend());
}
//Create output
TLorentzVector tmp; 
uint j=0; while (j<tmpoutput.size()){
   tmp.SetPtEtaPhiM( get<0>(tmpoutput[j]), get<1>(tmpoutput[j]), get<2>(tmpoutput[j]), 0.1056583745);
   output.push_back( make_tuple(tmp, get<3>(tmpoutput[j]),get<4>(tmpoutput[j]) )   );
   j++;
}
return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> SelectTracks(vector<float>* P1, vector<float>* P2, vector<float>* P3, vector<float>* chi2, vector<int>* nstubs, vector<int>* charge, vector<float>* z, bool isEMTF, bool isTT, bool overlapremoval,int endcap, CfgParser config)
{
// Note this function is designed to select stubs properties P4, CHARGE, z0
//Define variables
float mineta=stof(config.readStringOpt("L1TPreSelection::L1TMINETA"));
float maxeta=stof(config.readStringOpt("L1TPreSelection::L1TMAXETA"));
float minpt =stof(config.readStringOpt("L1TPreSelection::L1TMINPT"));
std::vector<std::tuple<float,float,float,int,float>> tmpoutput,tmpoutput2;
std::vector<std::tuple<TLorentzVector,int,float>> output;
//Make a tuple and remove overlapping track information 
std::vector<std::tuple<float,float,float,int,float>> info;
for(uint k=0;k<P3->size();k++)
{
  if(isEMTF) info.push_back(make_tuple(P1->at(k),P2->at(k),(P3->at(k)*TMath::Pi()/180), charge->at(k),0) );
  else if(isTT)
  {
  if(chi2->at(k) < 100 && nstubs->at(k)>=4) info.push_back(make_tuple(P1->at(k),P2->at(k),P3->at(k),charge->at(k),z->at(k) ));
  } 
  else
  {
   info.push_back(make_tuple(P1->at(k),P2->at(k),P3->at(k),charge->at(k),z->at(k) ));         
  }  
}
//Remove overlap of tracks
if(overlapremoval) info = RemoveOvelapppingTracks(info);
//If no tracks return empty
if(info.size()==0) return output;
//Fill inclusive endcap EMTF muons (for general studies)
if(endcap==0) 
{
    //Fill with positive endcap EMTF muons and ordered them by pt
    uint i=0; 
    while (i<info.size()) 
    {
      if(get<1>(info[i])>mineta && get<1>(info[i])<maxeta && get<0>(info[i])>minpt)
      {            
            tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
      }
      i++;
    }
    std::sort(tmpoutput.rbegin(), tmpoutput.rend());
    //Fill with negative endcap EMTF muons and ordered them by pt
    uint h=0; 
    while (h<info.size()) 
    {
      if(get<1>(info[h])<-mineta && get<1>(info[h])>-maxeta && get<0>(info[h])>minpt)
      {            
            tmpoutput2.push_back(make_tuple(get<0>(info[h]),get<1>(info[h]),get<2>(info[h]),get<3>(info[h]),get<4>(info[h])  ));
      }
      h++;
    }
    std::sort(tmpoutput2.rbegin(), tmpoutput2.rend());
    //Mergue both endcaps
    tmpoutput.insert(tmpoutput.end(),tmpoutput2.begin(),tmpoutput2.end());
}
else
{
    uint i=0; 
    while (i<info.size()) 
    {
      if(endcap>0) //Positive endcap (for signal only)
      {
         if(get<1>(info[i])>mineta && get<1>(info[i])<maxeta && get<0>(info[i])>minpt)
         {            
               tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
         }
         i++;
      }
      else        //Negative endcap (for signal only)
      {
         if(get<1>(info[i])<-mineta && get<1>(info[i])>-maxeta && get<0>(info[i])>minpt)
         {            
               tmpoutput.push_back(make_tuple(get<0>(info[i]),get<1>(info[i]),get<2>(info[i]),get<3>(info[i]),get<4>(info[i])  ));
         }
         i++;
      }
    }
    //Order by pt
    std::sort(tmpoutput.rbegin(), tmpoutput.rend());
}
//Create output
TLorentzVector tmp; 
uint j=0; while (j<tmpoutput.size()){
   tmp.SetPtEtaPhiM( get<0>(tmpoutput[j]), get<1>(tmpoutput[j]), get<2>(tmpoutput[j]), 0.1056583745);
   output.push_back( make_tuple(tmp, get<3>(tmpoutput[j]),get<4>(tmpoutput[j]) )   );
   j++;
}
return output;
}

std::vector<TLorentzVector> BuildVectorOfTLorentzVector(vector<float>* P1, vector<float>* P2, vector<float>* P3,vector<float>* P4)
{

//Order by PT
std::vector<TLorentzVector> output;

if(P1->size()==0) return output;

std::vector<std::tuple<float,float,float,float>> info;
uint j=0; while (j<P1->size())
{
   info.push_back(make_tuple(P1->at(j),P2->at(j),P3->at(j),P4->at(j)));
   j++;
}
std::sort(info.rbegin(), info.rend());
//Create output
TLorentzVector tmp; 
uint k=0; while (k<P1->size())
{
   tmp.SetPtEtaPhiE(get<0>(info[k]), get<1>(info[k]), get<2>(info[k]), get<3>(info[k]));
   output.push_back(tmp);
   k++;
}

return output;
}

std::vector<std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>>> SelectParticles(std::vector<TLorentzVector> taus, std::vector<TLorentzVector> muons, CfgParser config)
{
//Parse selection cuts from config file
uint NTAUS=stoi(config.readStringOpt("Generator::NTAUS"));
float TAUMINP  =stof(config.readStringOpt("Generator::TAUMINP"));
float TAUMINPT =stof(config.readStringOpt("Generator::TAUMINPT"));
uint  NMUS=stoi(config.readStringOpt("Generator::NMUS"));
float  MUMINP  =stof(config.readStringOpt("Generator::MUMINP"));
float  MUMINPT =stof(config.readStringOpt("Generator::MUMINPT"));
float  MUMINETA=stof(config.readStringOpt("Generator::MUMINETA"));
float  MUMAXETA=stof(config.readStringOpt("Generator::MUMAXETA"));
uint  NFIDMUS=stoi(config.readStringOpt("Generator::NFIDMUS"));
float  MUFIDMINPT=stof(config.readStringOpt("Generator::MUFIDMINPT"));
float  MUFIDMINETA=stof(config.readStringOpt("Generator::MUFIDMINETA"));
float  MUFIDMAXETA=stof(config.readStringOpt("Generator::MUFIDMAXETA"));

std::vector<std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>>> output;
//Check number of taus is one
if(taus.size()    != NTAUS) return output;
if(taus.at(0).P() <= TAUMINP || taus.at(0).Pt() <= TAUMINPT ) return output;
//Count muons satisfying the MC conditions
uint poscut1=0,negcut1=0;
for(uint k=0;k<3;k++)
{
   if(muons.at(k).Eta() > MUMINETA     && muons.at(k).Eta()<MUMAXETA      && muons.at(k).Pt()>MUMINPT) poscut1++;
   if(muons.at(k).Eta() <(-1*MUMINETA) && muons.at(k).Eta()>(-1*MUMAXETA) && muons.at(k).Pt()>MUMINPT) negcut1++;
}
if(poscut1<NMUS && negcut1<NMUS) return output;
uint poscut2=0,negcut2=0;
for(uint k=0;k<3;k++)
{
   if(muons.at(k).Eta() > MUMINETA     && muons.at(k).Eta()< MUMAXETA      && muons.at(k).P()> MUMINP)  poscut2++;
   if(muons.at(k).Eta() <(-1*MUMINETA) && muons.at(k).Eta()> (-1*MUMAXETA) && muons.at(k).P()> MUMINP)  negcut2++;
}
if(poscut2<NMUS && negcut2<NMUS) return output;
//-----SPACE FOR ADDITIONAL CUTS
uint poscut3=0,negcut3=0;
for(uint k=0;k<3;k++)
{
   if(muons.at(k).Eta() > MUFIDMINETA     && muons.at(k).Eta()   < MUFIDMAXETA  && muons.at(k).Pt() > MUFIDMINPT) poscut3++;
   if(muons.at(k).Eta() <(-1*MUFIDMINETA) && muons.at(k).Eta()>(-1*MUFIDMAXETA) && muons.at(k).Pt() > MUFIDMINPT) negcut3++;
}
if(poscut3<NFIDMUS && negcut3<NFIDMUS) return output;
//Create output
output.push_back(make_tuple(taus,muons));
//Return output
return output;
}

std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> GetObjectTracks(std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::array<float,11>>>> objects, int objectid)
{
 std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> object;
 uint j=0; while (j<objects.size())
 {
  if (objectid==0) object.push_back(get<2>(objects[j]));
  if (objectid==1) object.push_back(get<3>(objects[j])); 
  if (objectid==2) object.push_back(get<4>(objects[j]));
  if (objectid==3) object.push_back(get<5>(objects[j])); 
  if (objectid==4) object.push_back(get<6>(objects[j])); 
  j++;
 } 
 return object;
}

std::vector<std::vector<TLorentzVector>> GetObjectGens(std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::array<float,11>>>> objects, int objectid)
{
 std::vector<std::vector<TLorentzVector>> object;
 uint j=0; while (j<objects.size())
 {
  if (objectid==0) object.push_back(get<0>(objects[j]));
  if (objectid==1) object.push_back(get<1>(objects[j]));  
  j++;
 } 
 return object;
}

std::vector<std::vector<std::array<float,11>>> GetObjectHits(std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::array<float,11>> >> objects)
{
 std::vector<std::vector<std::array<float,11>>> object;
 uint j=0; while (j<objects.size())
 {
  object.push_back(get<7>(objects[j])); 
  j++;
 } 
 return object;
}

std::vector<std::string> split(std::string str,std::string sep){
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current=strtok(cstr,sep.c_str());
    while(current!=NULL){
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;
}


string createxlabel(string histoname){
   std::vector<std::string> arraylabel =split(Form("%s",histoname.c_str()),"_");

   string label="";
   for (uint k=0;k<arraylabel.size();k++)
   {
   //Importance
   if(arraylabel[k]=="0") label+="1st "; 
   if(arraylabel[k]=="1") label+="2nd ";   
   if(arraylabel[k]=="2") label+="3rd "; 
   }

   for (uint k=0;k<arraylabel.size();k++)
   {
   //Level
   cout<<arraylabel[k]<<endl;
   if(arraylabel[k]=="gen") label+="Gen-"; 
   if(arraylabel[k]=="rec") label+="Reco-";  
   //Particle
   if(arraylabel[k]=="mu") label+="#mu ";    
   if(arraylabel[k]=="tau") label+="#tau "; 
   //Kinematic variable
   if(arraylabel[k]=="p") label+="P [GeV]"; 
   if(arraylabel[k]=="pt") label+="P_{T} [GeV]";    
   if(arraylabel[k]=="eta") label+="#eta";
   }

   return label;
}