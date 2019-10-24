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
#include <string>
#include "CfgParser.h"

using namespace std;

std::vector<std::tuple<TLorentzVector,int,float>> FindMixedTrackCandidates(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12,std::vector<std::tuple<TLorentzVector,int,float>> Tracks3)
{
   std::vector<std::tuple<TLorentzVector,int,float>> output,postracks,negtracks;
   //First look for tracks in the same endcap
   for(uint k=0; k<Tracks12.size(); k++)
   { 
     if( get<0>(Tracks12.at(k)).Eta() > 0) postracks.push_back(Tracks12.at(k));
     if( get<0>(Tracks12.at(k)).Eta() < 0) negtracks.push_back(Tracks12.at(k));
   }

   //If found at least two tracks, then look for tracks compatible with the tau charge
   //Check candidate in positive endcap
   if(postracks.size()>=2)
   {
         for(uint i=0;i<2;i++) output.push_back(postracks.at(i) );
         for(uint k=0; k<Tracks3.size(); k++)
         {
            int charge = get<1>(output.at(0)) + get<1>(output.at(1)) + get<1>(Tracks3.at(k));
            if(get<0>(Tracks3.at(k)).Eta() > 0 && abs(charge)==1 ) output.push_back(Tracks3.at(k));
         } 
   }
   else
   {
     //Check candidate in negative endcap
     if(negtracks.size()>=2)
     {
           for(uint i=0;i<2;i++) output.push_back(negtracks.at(i));
           for(uint k=0; k<Tracks3.size(); k++)
           {
              int charge = get<1>(output.at(0)) + get<1>(output.at(1)) + get<1>(Tracks3.at(k));
              if( get<0>(Tracks3.at(k)).Eta() < 0 && abs(charge)==1) output.push_back(Tracks3.at(k));
           }
     }
   }

   return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> FindPureTrackCandidates(std::vector<std::tuple<TLorentzVector,int,float>> Tracks)
{
   std::vector<std::tuple<TLorentzVector,int,float>> output,postracks,negtracks;
   bool triggeredpos=false,triggeredneg=false;
   //First get the two leading pt tracks of the set
   int pos2=0,neg2=0,kpos=0,kneg=0;
   for(uint k=0; k<Tracks.size(); k++)
   { 
     if( get<0>(Tracks.at(k)).Eta() > 0 && pos2<2 ){postracks.push_back(Tracks.at(k));pos2++;kpos=k;}
     if( get<0>(Tracks.at(k)).Eta() < 0 && neg2<2 ){negtracks.push_back(Tracks.at(k));neg2++;kneg=k;}
   }
   //Now just fill with the one that have the right charge and trigger it
   if(pos2==2)
     for(uint k=kpos+1; k<Tracks.size(); k++)
     { 
       int charge = get<1>(postracks.at(0)) + get<1>(postracks.at(1)) + get<1>(Tracks.at(k));
       if( get<0>(Tracks.at(k)).Eta() > 0 && abs(charge)==1 ){postracks.push_back(Tracks.at(k));triggeredpos=true;}
     }
   if(neg2==2)
     for(uint k=kneg+1; k<Tracks.size(); k++)
     { 
       int charge = get<1>(negtracks.at(0)) + get<1>(negtracks.at(1)) + get<1>(Tracks.at(k));
       if( get<0>(Tracks.at(k)).Eta() < 0 && abs(charge)==1 ){negtracks.push_back(Tracks.at(k));triggeredneg=true;}
     }
   //Save information
   if(       triggeredpos && !triggeredneg ) for(uint i=0;i<postracks.size();i++) output.push_back(postracks.at(i)); 
   else if (!triggeredpos && triggeredneg  ) for(uint i=0;i<negtracks.size();i++) output.push_back(negtracks.at(i));
   else if ( triggeredpos && triggeredneg  ) for(uint i=0;i<postracks.size();i++) output.push_back(postracks.at(i));
   else
   {
   //Do nothing
   }

   return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> CompleteInvariantMassDistanceTrack(std::vector<std::tuple<TLorentzVector,int,float>> tracks, float mwindow, float drwindow, float dzwindow)
{
bool found=false;
std::vector<std::tuple<TLorentzVector,int,float>> output;
std::tuple<TLorentzVector,int,float> track3;
//Check if there is a track that minimizes the tau invariant mass and window  
float range=mwindow,dM;
for(uint i=2; i<tracks.size(); ++i)
{
  dM = abs(( get<0>(tracks.at(0)) + get<0>(tracks.at(1)) + get<0>(tracks.at(i)) ).M()-1.77); 
  if(dM<range){found=true;track3=tracks.at(i); range=dM;}
}
//Found the candidate?
if(!found) return output;
//The candidate is within the dR window?
float maxdR;
if(drwindow != 0)
  maxdR = std::max( {get<0>(tracks.at(0)).DeltaR(get<0>(tracks.at(1))), get<0>(tracks.at(0)).DeltaR(get<0>(track3)), get<0>(tracks.at(1)).DeltaR(get<0>(track3))});
else
  maxdR = -1;
if(maxdR>=drwindow) return output;
//The candidate is within the dZ window?
float maxdz;
if(dzwindow != 0)
  maxdz = std::max( { abs(get<2>(tracks.at(0))-get<2>(tracks.at(1))), abs(get<2>(tracks.at(0))-get<2>(track3) ), abs(get<2>(tracks.at(1))-get<2>(track3)) });
else
  maxdz = -1;
if(maxdz>=dzwindow) return output;
//Build output
output.push_back( tracks.at(0) );
output.push_back( tracks.at(1) );
output.push_back( track3       );

return output;
}


std::vector<std::tuple<TLorentzVector,int,float>> CompleteME0Track(std::vector<std::tuple<TLorentzVector,int,float>> tracks, float dzwindow, float drwindow, float bendcut)
{
//Define output 
std::vector<std::tuple<TLorentzVector,int,float>> output;
//First cut on dz between the selected two tracker-based tracks
float maxdz;
if(dzwindow != 0)
  maxdz = abs(get<2>(tracks.at(0))-get<2>(tracks.at(1)));
else
  maxdz = -1;
if(maxdz>=dzwindow) return output; 
//Select track that minized dR 
bool found=false;
std::tuple<TLorentzVector,int,float> track3;
//Check if there is a track candidate
float range=drwindow,maxdR;
for(uint i=2; i<tracks.size(); ++i)
{ 
  //Cut on eta w.r.t. to the tracks
  if( abs(get<0>(tracks.at(i)).Eta()) <= abs(get<0>(tracks.at(0)).Eta()) ) continue;
  if( abs(get<0>(tracks.at(i)).Eta()) <= abs(get<0>(tracks.at(1)).Eta()) ) continue; 
  //Cut on eta for ME0 stubs
  //if( abs(get<0>(tracks.at(i)).Eta()) < 2.0) continue;
  //Check for minimizing dR
  maxdR = std::max( {  get<0>(tracks.at(0)).DeltaR(get<0>(tracks.at(1))), get<0>(tracks.at(0)).DeltaR(get<0>(tracks.at(i))), get<0>(tracks.at(1)).DeltaR(get<0>(tracks.at(i)))});
  if(maxdR < range ){found=true;track3=tracks.at(i); range=maxdR;}
}
//Found the candidate?
if(!found) return output;
//Build output
output.push_back( tracks.at(0) );
output.push_back( tracks.at(1) );
output.push_back( track3       );

return output;
}



std::vector<std::tuple<TLorentzVector,int,float>> CompleteDistanceTrack(std::vector<std::tuple<TLorentzVector,int,float>> tracks, float drwindow, float dzwindow)
{
bool found=false;
std::vector<std::tuple<TLorentzVector,int,float>> output;
std::tuple<TLorentzVector,int,float> track3;
//Check if there is a track that minimizes the tau invariant mass and window  
float range=drwindow,maxdR;
for(uint i=2; i<tracks.size(); ++i)
{ 
  maxdR = std::max( {  get<0>(tracks.at(0)).DeltaR(get<0>(tracks.at(1))), get<0>(tracks.at(0)).DeltaR(get<0>(tracks.at(i))), get<0>(tracks.at(1)).DeltaR(get<0>(tracks.at(i)))});
  if(maxdR < range){found=true;track3=tracks.at(i); range=maxdR;}
}
//Found the candidate?
if(!found) return output;
//The candidate is within the dZ window?
float maxdz;
if(dzwindow != 0)
  maxdz = std::max( { abs(get<2>(tracks.at(0))-get<2>(tracks.at(1))), abs(get<2>(tracks.at(0))-get<2>(track3) ), abs(get<2>(tracks.at(1))-get<2>(track3)) });
else
  maxdz = -1;
if(maxdz>=dzwindow) return output;
//Is candidate within a reasonable invariant mass?
float dM = abs((  get<0>(tracks.at(0))  + get<0>( tracks.at(1) ) + get<0>(track3)  ).M()-1.77); 
if(dM >= 0.3) return output;
//Build output
output.push_back(tracks.at(0));
output.push_back(tracks.at(1));
output.push_back(track3      );

return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> CompleteTrack(std::vector<std::tuple<TLorentzVector,int,float>> tracks)
{
std::vector<std::tuple<TLorentzVector,int,float>> output;
//Build output
output.push_back(tracks.at(0) );
output.push_back(tracks.at(1) );
output.push_back(tracks.at(2) );

return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> MixedME0TriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12, std::vector<std::tuple<TLorentzVector,int,float>> Tracks3, float dzwindow, float drwindow, float bendcut)
{
    std::vector<std::tuple<TLorentzVector,int,float>> output,tracks;
    //Check for candiDates in positive or negative endcap
    tracks=FindMixedTrackCandidates(Tracks12,Tracks3);
    //Check that there are three tracks
    if(tracks.size() < 3) return output;
    //Select the three tracks based on the configuration
    output = CompleteME0Track(tracks, dzwindow, drwindow, bendcut);
    //return output
    return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> MixedTriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12, std::vector<std::tuple<TLorentzVector,int,float>> Tracks3,float mwindow, float drwindow, float dzwindow)
{
    std::vector<std::tuple<TLorentzVector,int,float>> output,tracks;
    //Check for candiDates in positive or negative endcap
    tracks=FindMixedTrackCandidates(Tracks12,Tracks3);
    //Check that there are three tracks
    if(tracks.size() < 3) return output;
    //Select the three tracks based on the configuration
    if(mwindow != 0)
      output = CompleteInvariantMassDistanceTrack(tracks, mwindow, drwindow, dzwindow);
    else if (mwindow ==0 && drwindow != 0)
      output = CompleteDistanceTrack(tracks, drwindow, dzwindow);
    else 
      output = CompleteTrack(tracks);
    //return output
    return output;
}

std::vector<std::tuple<TLorentzVector,int,float>> PureTriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks,float mwindow, float drwindow, float dzwindow)
{
    std::vector<std::tuple<TLorentzVector,int,float>> output,tracks;
    //Check for candiDates in positive or negative endcap
    tracks=FindPureTrackCandidates(Tracks);
    //Check that there are three tracks
    if(tracks.size() < 3) return output;
    //Select the three tracks based on the configuration
    if(mwindow != 0)
      output = CompleteInvariantMassDistanceTrack(tracks, mwindow, drwindow, dzwindow);
    else if (mwindow ==0 && drwindow != 0)
      output = CompleteDistanceTrack(tracks, drwindow, dzwindow);
    else 
      output = CompleteTrack(tracks);
    //return output
    return output;
}

void MixedTripleME0TriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12, std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float dzwindow, float drwindow, float bendcut, string tagname)
{
//Book Histogram for the reference
TH1F *h_TRG_MIXED_maxdR        = new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "",100, 0,  4);
TH1F *h_TRG_MIXED_totalcharge  = new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG_MIXED_mu_dz0_01    = new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_MIXED_bend         = new TH1F(Form("h_TRG%s_bend",  tagname.c_str() ),  "", 50,  0, 0.5);
TH1F *h_TRG_MIXED_inversebend  = new TH1F(Form("h_TRG%s_inversebend",  tagname.c_str() ),  "", 100,  0, 100);
TH1F *h_TRG_MIXED_variable     = new TH1F(Form("h_TRG_%s_variable",tagname.c_str()),"" ,       400,  0, 100);
//Fill variable histogram
int trigger=0;
//int triggerplus=0;
for (uint i=0; i<Tracks12.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = MixedME0TriggerSelector(Tracks12.at(i),Tracks3.at(i), dzwindow, drwindow, bendcut);
  if(triggered.size() == 3)
  {    
    h_TRG_MIXED_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
    h_TRG_MIXED_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
    h_TRG_MIXED_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
    h_TRG_MIXED_bend->Fill(  abs(get<2>(triggered.at(2)) )  );
    h_TRG_MIXED_inversebend->Fill( 1 / abs(get<2>(triggered.at(2)) ) );
    h_TRG_MIXED_variable->Fill(    1 / abs(get<2>(triggered.at(2)) ) );
    trigger++;
  } 
  else
  {
  h_TRG_MIXED_variable->Fill(-1);
  }
}

MakeRatePlot(h_TRG_MIXED_variable, 2760*11.246,"rate");
MakeRatePlot(h_TRG_MIXED_variable, 1,"roc");
//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Acceptance for Trigger %s   = ",tagname.c_str() )<<(float)trigger/Tracks12.size()<<endl;
cout<<Form("Maximum Rate for Trigger %s   (kHz) = ",tagname.c_str() )<<(float)2760*11.246*trigger/Tracks12.size()<<endl;
}

void MixedTripleMuonTriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12, std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float mwindow, float drwindow, float dzwindow, string tagname)
{
//Book Histogram for the reference
TH1F *h_TRG_MIXED_variable = new TH1F(Form("h_TRG_%s_variable",tagname.c_str()),"" ,304, 2, 40);
TH1F *h_TRG_MIXED_invmass  = new TH1F(Form("h_TRG%s_invmass",  tagname.c_str() ),"",100, 0,  3);
TH1F *h_TRG_MIXED_maxdR    = new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "",100, 0,  4);
TH1F *h_TRG_MIXED_totalcharge    =  new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG_MIXED_mu_dz0_01      =  new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_MIXED_mu_dz0_12      =  new TH1F(Form("h_TRG%s_mu_dz0_12",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_MIXED_mu_dz0_02      =  new TH1F(Form("h_TRG%s_mu_dz0_02",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_MIXED_mu_dz0_max     =  new TH1F(Form("h_TRG%s_mu_dz0_max",tagname.c_str() ), "", 80, 0, 20);
//Fill variable histogram
int trigger=0;
//int triggerplus=0;
for (uint i=0; i<Tracks12.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = MixedTriggerSelector(Tracks12.at(i),Tracks3.at(i), mwindow, drwindow, dzwindow);
  if(triggered.size() == 3)
  {    
    h_TRG_MIXED_variable->Fill(get<0>(triggered.at(1)).Pt());
    h_TRG_MIXED_invmass->Fill( (get<0>(triggered.at(0))+get<0>(triggered.at(1))+get<0>(triggered.at(2))).M() );
    h_TRG_MIXED_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
    h_TRG_MIXED_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
    h_TRG_MIXED_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
    h_TRG_MIXED_mu_dz0_12->Fill( abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) )  );
    h_TRG_MIXED_mu_dz0_02->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )  );
    h_TRG_MIXED_mu_dz0_max->Fill( std::max( { abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) ) , abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) ) , abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )})    );       
    trigger++;
  } 
  else
  {
  h_TRG_MIXED_variable->Fill(-1);
  }
}

MakeRatePlot(h_TRG_MIXED_variable, 2760*11.246,"rate");
MakeRatePlot(h_TRG_MIXED_variable, 1,"roc");
//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Acceptance for Trigger %s   = ",tagname.c_str() )<<(float)trigger/Tracks12.size()<<endl;
cout<<Form("Maximum Rate for Trigger %s   (kHz) = ",tagname.c_str() )<<(float)2760*11.246*trigger/Tracks12.size()<<endl;
}

void PureTripleMuonTriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks,  float mwindow, float drwindow, float dzwindow, string tagname)
{
//Book Histogram for the reference
TH1F *h_TRG_PURE_variable       =  new TH1F(Form("h_TRG_%s_variable",tagname.c_str()), "", 304, 2, 40);
TH1F *h_TRG_PURE_invmass        =  new TH1F(Form("h_TRG%s_invmass",  tagname.c_str() ),"", 100, 0,  3);
TH1F *h_TRG_PURE_maxdR          =  new TH1F(Form("h_TRG%s_maxdR",  tagname.c_str() ),  "", 100, 0,  4);
TH1F *h_TRG_PURE_totalcharge    =  new TH1F(Form("h_TRG%s_totalcharge",tagname.c_str() ), "", 8, -4,4);
TH1F *h_TRG_PURE_mu_dz0_01      =  new TH1F(Form("h_TRG%s_mu_dz0_01",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_PURE_mu_dz0_12      =  new TH1F(Form("h_TRG%s_mu_dz0_12",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_PURE_mu_dz0_02      =  new TH1F(Form("h_TRG%s_mu_dz0_02",tagname.c_str() ), "", 80, 0, 20);
TH1F *h_TRG_PURE_mu_dz0_max     =  new TH1F(Form("h_TRG%s_mu_dz0_max",tagname.c_str() ), "", 80, 0, 20);
//Fill variable histogram
int trigger=0;
//int triggerplus=0;
for (uint i=0; i<Tracks.size(); ++i)
{  
  std::vector<std::tuple<TLorentzVector,int,float>> triggered = PureTriggerSelector(Tracks.at(i), mwindow, drwindow, dzwindow);
  if(triggered.size() == 3)
  {
    h_TRG_PURE_variable->Fill(get<0>(triggered.at(1)).Pt());    
    h_TRG_PURE_invmass->Fill( (get<0>(triggered.at(0))+get<0>(triggered.at(1))+get<0>(triggered.at(2))).M() );
    h_TRG_PURE_maxdR->Fill( std::max( {get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(1))), get<0>(triggered.at(0)).DeltaR(get<0>(triggered.at(2))) , get<0>(triggered.at(1)).DeltaR(get<0>(triggered.at(2))  )})   );
    h_TRG_PURE_totalcharge->Fill( get<1>(triggered.at(0)) + get<1>(triggered.at(1)) + get<1>(triggered.at(2)) );
    h_TRG_PURE_mu_dz0_01->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) )  );
    h_TRG_PURE_mu_dz0_12->Fill( abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) )  );
    h_TRG_PURE_mu_dz0_02->Fill( abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )  ); 
    h_TRG_PURE_mu_dz0_max->Fill( std::max( { abs( get<2>(triggered.at(0)) - get<2>(triggered.at(1)) ) , abs( get<2>(triggered.at(1)) - get<2>(triggered.at(2)) ) , abs( get<2>(triggered.at(0)) - get<2>(triggered.at(2)) )})    ); 
    trigger++;
  } 
  else
  {
  h_TRG_PURE_variable->Fill(-1);
  }
}

//Make rate plot
MakeRatePlot(h_TRG_PURE_variable, 2760*11.246,"rate");
MakeRatePlot(h_TRG_PURE_variable, 1,"roc");
//REPORT
cout<<"TRIGGER REPORT:"<<endl;
cout<<Form("Trigger Acceptance for Trigger %s   = ",tagname.c_str() )<<(float)trigger/Tracks.size()<<endl;
cout<<Form("Maximum Rate for Trigger %s   (kHz) = ",tagname.c_str() )<<(float)2760*11.246*trigger/Tracks.size()<<endl;
}


//Matching muons to track modules
void MatchingCalculator(std::vector<std::vector<TLorentzVector>> Muons,std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks, string tagname)
{
  //Define histograms here
  //MAKE BINS
  Float_t etabins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.5,4.0}; 
  int netabins = sizeof(etabins) / sizeof(etabins[0])-1;
  Float_t ptbins[] = {0,0.5,1.0,1.5,2,2.5,3,4,5,15};
  int nptbins = sizeof(ptbins) / sizeof(ptbins[0])-1;
  Float_t pbins[] = {0,2.5,5,10,15,20,30,100};
  int npbins = sizeof(pbins) / sizeof(pbins[0])-1;
  TH1F *h_MCH_gen_mu_p_0     =  new TH1F(Form("h_MCH_%s_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
  TH1F *h_MCH_gen_mu_pt_0    =  new TH1F(Form("h_MCH_%s_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_MCH_gen_mu_eta_0   =  new TH1F(Form("h_MCH_%s_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
  TH1F *h_MCH_gen_mu_p_1     =  new TH1F(Form("h_MCH_%s_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
  TH1F *h_MCH_gen_mu_pt_1    =  new TH1F(Form("h_MCH_%s_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_MCH_gen_mu_eta_1   =  new TH1F(Form("h_MCH_%s_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
  TH1F *h_MCH_gen_mu_p_2     =  new TH1F(Form("h_MCH_%s_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
  TH1F *h_MCH_gen_mu_pt_2    =  new TH1F(Form("h_MCH_%s_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_MCH_gen_mu_eta_2   =  new TH1F(Form("h_MCH_%s_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins); 
  TH2F *h_MCH_gen_mu_pt_eta_0 =  new TH2F(Form("h_MCH_%s_gen_mu_pt_eta_0",tagname.c_str()), "", npbins,pbins,netabins,etabins);
  TH2F *h_MCH_gen_mu_pt_eta_1 =  new TH2F(Form("h_MCH_%s_gen_mu_pt_eta_1",tagname.c_str()), "", npbins,pbins,netabins,etabins);
  TH2F *h_MCH_gen_mu_pt_eta_2 =  new TH2F(Form("h_MCH_%s_gen_mu_pt_eta_2",tagname.c_str()), "", 20,0,20,20,1,3);
  TH2F *h_MCH_gen_mu_n_eta_2  =  new TH2F(Form("h_MCH_%s_gen_mu_n_eta_2",tagname.c_str()), "", 10,0,10,20,1,3);
  TH2F *h_MCH_gen_mu_pt_n_2   =  new TH2F(Form("h_MCH_%s_gen_mu_pt_n_2",tagname.c_str()), "",  20,0,20,10,0,10);

  TH1F *h_NOMCH_gen_mu_p_0   =  new TH1F(Form("h_NOMCH_%s_gen_mu_p_0",tagname.c_str()), "", npbins,pbins);
  TH1F *h_NOMCH_gen_mu_pt_0  =  new TH1F(Form("h_NOMCH_%s_gen_mu_pt_0",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_NOMCH_gen_mu_eta_0 =  new TH1F(Form("h_NOMCH_%s_gen_mu_eta_0",tagname.c_str()),"", netabins, etabins);
  TH1F *h_NOMCH_gen_mu_p_1   =  new TH1F(Form("h_NOMCH_%s_gen_mu_p_1",tagname.c_str()), "", npbins,pbins);
  TH1F *h_NOMCH_gen_mu_pt_1  =  new TH1F(Form("h_NOMCH_%s_gen_mu_pt_1",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_NOMCH_gen_mu_eta_1 =  new TH1F(Form("h_NOMCH_%s_gen_mu_eta_1",tagname.c_str()),"", netabins, etabins);
  TH1F *h_NOMCH_gen_mu_p_2   =  new TH1F(Form("h_NOMCH_%s_gen_mu_p_2",tagname.c_str()), "", npbins,pbins);
  TH1F *h_NOMCH_gen_mu_pt_2  =  new TH1F(Form("h_NOMCH_%s_gen_mu_pt_2",tagname.c_str()), "",nptbins,ptbins);
  TH1F *h_NOMCH_gen_mu_eta_2 =  new TH1F(Form("h_NOMCH_%s_gen_mu_eta_2",tagname.c_str()),"", netabins, etabins); 
  TH2F *h_NOMCH_gen_mu_pt_eta_0 =  new TH2F(Form("h_NOMCH_%s_gen_mu_pt_eta_0",tagname.c_str()), "", npbins,pbins,netabins,etabins);
  TH2F *h_NOMCH_gen_mu_pt_eta_1 =  new TH2F(Form("h_NOMCH_%s_gen_mu_pt_eta_1",tagname.c_str()), "", npbins,pbins,netabins,etabins);
  TH2F *h_NOMCH_gen_mu_pt_eta_2 =  new TH2F(Form("h_NOMCH_%s_gen_mu_pt_eta_2",tagname.c_str()), "", 20,0,20,20,1,3);
  TH2F *h_NOMCH_gen_mu_n_eta_2 =  new TH2F(Form("h_NOMCH_%s_gen_mu_n_eta_2",tagname.c_str()), "",   10,0,10,20,1,3);
  TH2F *h_NOMCH_gen_mu_pt_n_2  =  new TH2F(Form("h_NOMCH_%s_gen_mu_pt_n_2",tagname.c_str()), "",    20,0,20,10,0,10);

  TH1F *h_MCH_gen_mu_ptres_0 =  new TH1F(Form("h_MCH_%s_gen_mu_ptres_0",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_gen_mu_ptres_1 =  new TH1F(Form("h_MCH_%s_gen_mu_ptres_1",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_gen_mu_ptres_2 =  new TH1F(Form("h_MCH_%s_gen_mu_ptres_2",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_gen_mu_dR_0    =  new TH1F(Form("h_MCH_%s_gen_mu_dR_0",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_gen_mu_dR_1    =  new TH1F(Form("h_MCH_%s_gen_mu_dR_1",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_gen_mu_dR_2    =  new TH1F(Form("h_MCH_%s_gen_mu_dR_2",tagname.c_str()), "",40,-2,2);
  TH1F *h_MCH_nmatches       =  new TH1F(Form("h_MCH_%s_nmatches"  ,tagname.c_str()), "", 4,0,4);  
  TH1F *h_MCH_invmass        =  new TH1F(Form("h_MCH_%s_invmass"  ,tagname.c_str()), "", 100,0,4); 


  TH1F *h_MCH_gen_mu_dR01_2   = new TH1F(Form("h_MCH_%s_gen_mu_dR01_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dR02_2   = new TH1F(Form("h_MCH_%s_gen_mu_dR02_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dR12_2   = new TH1F(Form("h_MCH_%s_gen_mu_dR12_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dRmin_2  = new TH1F(Form("h_MCH_%s_gen_mu_dRmin_2" ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dEta01_2 = new TH1F(Form("h_MCH_%s_gen_mu_dEta01_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dEta02_2 = new TH1F(Form("h_MCH_%s_gen_mu_dEta02_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dEta12_2 = new TH1F(Form("h_MCH_%s_gen_mu_dEta12_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_MCH_gen_mu_dEtamin_2  = new TH1F(Form("h_MCH_%s_gen_mu_dEtamin_2",tagname.c_str()),"", 40, 0, 1);
  TH1F *h_MCH_gen_mu_dPhi01_2   = new TH1F(Form("h_MCH_%s_gen_mu_dPhi01_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_MCH_gen_mu_dPhi02_2   = new TH1F(Form("h_MCH_%s_gen_mu_dPhi02_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_MCH_gen_mu_dPhi12_2   = new TH1F(Form("h_MCH_%s_gen_mu_dPhi12_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_MCH_gen_mu_dPhimin_2  = new TH1F(Form("h_MCH_%s_gen_mu_dPhimin_2",tagname.c_str()),"", 40, 0, 1);
  TH1F *h_MCH_ntracks_0  = new TH1F(Form("h_MCH_%s_ntracks_0",tagname.c_str()),"", 10, 0, 10);
  TH1F *h_MCH_ntracks_1  = new TH1F(Form("h_MCH_%s_ntracks_1",tagname.c_str()),"", 10, 0, 10);
  TH1F *h_MCH_ntracks_2  = new TH1F(Form("h_MCH_%s_ntracks_2",tagname.c_str()),"", 10, 0, 10);
  TH1F *h_NOMCH_gen_mu_dR01_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dR01_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dR02_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dR02_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dR12_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dR12_2"  ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dRmin_2  = new TH1F(Form("h_NOMCH_%s_gen_mu_dRmin_2" ,tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dEta01_2 = new TH1F(Form("h_NOMCH_%s_gen_mu_dEta01_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dEta02_2 = new TH1F(Form("h_NOMCH_%s_gen_mu_dEta02_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dEta12_2 = new TH1F(Form("h_NOMCH_%s_gen_mu_dEta12_2",tagname.c_str() ), "",  40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dEtamin_2  = new TH1F(Form("h_NOMCH_%s_gen_mu_dEtamin_2",tagname.c_str()),"", 40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dPhi01_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dPhi01_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dPhi02_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dPhi02_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dPhi12_2   = new TH1F(Form("h_NOMCH_%s_gen_mu_dPhi12_2" ,tagname.c_str()),"", 40, 0, 1);
  TH1F *h_NOMCH_gen_mu_dPhimin_2  = new TH1F(Form("h_NOMCH_%s_gen_mu_dPhimin_2",tagname.c_str()),"", 40, 0, 1);
  TH1F *h_NOMCH_ntracks_0  = new TH1F(Form("h_NOMCH_%s_ntracks_0",tagname.c_str()),"", 10, 0, 10);
  TH1F *h_NOMCH_ntracks_1  = new TH1F(Form("h_NOMCH_%s_ntracks_1",tagname.c_str()),"", 10, 0, 10);
  TH1F *h_NOMCH_ntracks_2  = new TH1F(Form("h_NOMCH_%s_ntracks_2",tagname.c_str()),"", 10, 0, 10);


  //Loop over events
  for(uint i=0; i<Muons.size(); ++i)
  {
    std::vector<std::tuple<bool,std::tuple<TLorentzVector,int,float>>> flags = RunMatching(Muons.at(i),Tracks.at(i)); 
    //Fill histograms 
    int matches=0; 
    if( get<0>(flags.at(0)) )
    {
      h_MCH_ntracks_0->Fill( Tracks.at(i).size() );
      h_MCH_gen_mu_p_0->Fill(   Muons.at(i).at(0).P());
      h_MCH_gen_mu_pt_0->Fill(  Muons.at(i).at(0).Pt());
      h_MCH_gen_mu_eta_0->Fill( abs(Muons.at(i).at(0).Eta())   );
      h_MCH_gen_mu_ptres_0->Fill( (Muons.at(i).at(0).Pt() - get<0>(get<1>(flags.at(0))).Pt())/Muons.at(i).at(0).Pt()   );
      h_MCH_gen_mu_dR_0->Fill( Muons.at(i).at(0).DeltaR(  get<0>(get<1>(flags.at(0)))  )   );
      h_MCH_gen_mu_pt_eta_0->Fill( Muons.at(i).at(0).Pt(), abs(Muons.at(i).at(0).Eta()  ) );
      matches++;
    }
    else{
      h_NOMCH_ntracks_0->Fill( Tracks.at(i).size() );
      h_NOMCH_gen_mu_p_0->Fill(   Muons.at(i).at(0).P());
      h_NOMCH_gen_mu_pt_0->Fill(  Muons.at(i).at(0).Pt());
      h_NOMCH_gen_mu_eta_0->Fill( abs(Muons.at(i).at(0).Eta()) );      
      h_NOMCH_gen_mu_pt_eta_0->Fill( Muons.at(i).at(0).Pt(), abs(Muons.at(i).at(0).Eta()  ) );
    } 
    if( get<0>(flags.at(1)) )
    {
      h_MCH_ntracks_1->Fill( Tracks.at(i).size() );
      h_MCH_gen_mu_p_1->Fill(   Muons.at(i).at(1).P());
      h_MCH_gen_mu_pt_1->Fill(  Muons.at(i).at(1).Pt());
      h_MCH_gen_mu_eta_1->Fill( abs(Muons.at(i).at(1).Eta() )  );
      h_MCH_gen_mu_ptres_1->Fill( (Muons.at(i).at(1).Pt() - get<0>(get<1>(flags.at(1)) ).Pt())/Muons.at(i).at(1).Pt()   );
      h_MCH_gen_mu_dR_1->Fill( Muons.at(i).at(1).DeltaR(  get<0>(get<1>(flags.at(1)) )  )   );
      h_MCH_gen_mu_pt_eta_1->Fill( Muons.at(i).at(1).Pt(), abs(Muons.at(i).at(1).Eta()  ) );
      matches++;
    } 
    else{
      h_NOMCH_ntracks_1->Fill( Tracks.at(i).size() );
      h_NOMCH_gen_mu_p_1->Fill(   Muons.at(i).at(1).P());
      h_NOMCH_gen_mu_pt_1->Fill(  Muons.at(i).at(1).Pt());
      h_NOMCH_gen_mu_eta_1->Fill(abs( Muons.at(i).at(1).Eta() )  );      
      h_NOMCH_gen_mu_pt_eta_1->Fill( Muons.at(i).at(1).Pt(), abs(Muons.at(i).at(1).Eta()  ) );
    }
    if( get<0>(flags.at(2)) )
    {
      h_MCH_ntracks_2->Fill( Tracks.at(i).size() );
      h_MCH_gen_mu_p_2->Fill(   Muons.at(i).at(2).P());
      h_MCH_gen_mu_pt_2->Fill(  Muons.at(i).at(2).Pt());
      h_MCH_gen_mu_eta_2->Fill( abs(Muons.at(i).at(2).Eta()  )  );
      h_MCH_gen_mu_ptres_2->Fill( (Muons.at(i).at(2).Pt() - get<0>(get<1>(flags.at(2))  ).Pt())/Muons.at(i).at(2).Pt()   );
      h_MCH_gen_mu_dR_2->Fill( Muons.at(i).at(2).DeltaR(  get<0>(get<1>(flags.at(2))  )  )   );
      h_MCH_gen_mu_dR01_2->Fill( Muons.at(i).at(0).DeltaR(Muons.at(i).at(1)));
      h_MCH_gen_mu_dR02_2->Fill( Muons.at(i).at(0).DeltaR(Muons.at(i).at(2)));
      h_MCH_gen_mu_dR12_2->Fill( Muons.at(i).at(1).DeltaR(Muons.at(i).at(2)));
      h_MCH_gen_mu_dRmin_2->Fill(std::min( {Muons.at(i).at(0).DeltaR(Muons.at(i).at(1)), Muons.at(i).at(0).DeltaR(Muons.at(i).at(2)) , Muons.at(i).at(1).DeltaR(Muons.at(i).at(2)  )})   );
      h_MCH_gen_mu_dPhi01_2->Fill( Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(1)  ));
      h_MCH_gen_mu_dPhi02_2->Fill( Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(2)  ));
      h_MCH_gen_mu_dPhi12_2->Fill( Muons.at(i).at(1).DeltaPhi(Muons.at(i).at(2)  ));
      h_MCH_gen_mu_dPhimin_2->Fill(std::min( {Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(1)), Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(2)), Muons.at(i).at(1).DeltaPhi(Muons.at(i).at(2))})   );
      h_MCH_gen_mu_dEta01_2->Fill( abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(1).Eta() ) );
      h_MCH_gen_mu_dEta02_2->Fill( abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(2).Eta() ) );
      h_MCH_gen_mu_dEta12_2->Fill( abs( Muons.at(i).at(1).Eta() - Muons.at(i).at(2).Eta() ) );
      h_MCH_gen_mu_dEtamin_2->Fill(std::min( {abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(1).Eta() ) , abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(2).Eta() ) , abs( Muons.at(i).at(1).Eta() - Muons.at(i).at(2).Eta() )   })   );
      h_MCH_gen_mu_pt_eta_2->Fill( Muons.at(i).at(2).Pt(), abs(Muons.at(i).at(2).Eta()  ) );
      h_MCH_gen_mu_n_eta_2->Fill( Tracks.at(i).size(), abs(Muons.at(i).at(2).Eta()  ) );
      h_MCH_gen_mu_pt_n_2->Fill( Muons.at(i).at(2).Pt(), Tracks.at(i).size() );
      matches++;
    }
    else{
      h_NOMCH_ntracks_2->Fill( Tracks.at(i).size() );
      h_NOMCH_gen_mu_p_2->Fill(   Muons.at(i).at(2).P());
      h_NOMCH_gen_mu_pt_2->Fill(  Muons.at(i).at(2).Pt());
      h_NOMCH_gen_mu_eta_2->Fill( abs(Muons.at(i).at(2).Eta() )  );      
      h_NOMCH_gen_mu_dR01_2->Fill( Muons.at(i).at(0).DeltaR(Muons.at(i).at(1)));
      h_NOMCH_gen_mu_dR02_2->Fill( Muons.at(i).at(0).DeltaR(Muons.at(i).at(2)));
      h_NOMCH_gen_mu_dR12_2->Fill( Muons.at(i).at(1).DeltaR(Muons.at(i).at(2)));
      h_NOMCH_gen_mu_dRmin_2->Fill(std::min( {Muons.at(i).at(0).DeltaR(Muons.at(i).at(1)), Muons.at(i).at(0).DeltaR(Muons.at(i).at(2)) , Muons.at(i).at(1).DeltaR(Muons.at(i).at(2)  )})   );
      h_NOMCH_gen_mu_dPhi01_2->Fill( Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(1)  ));
      h_NOMCH_gen_mu_dPhi02_2->Fill( Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(2)  ));
      h_NOMCH_gen_mu_dPhi12_2->Fill( Muons.at(i).at(1).DeltaPhi(Muons.at(i).at(2)  ));
      h_NOMCH_gen_mu_dPhimin_2->Fill(std::min( {Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(1)), Muons.at(i).at(0).DeltaPhi(Muons.at(i).at(2)), Muons.at(i).at(1).DeltaPhi(Muons.at(i).at(2))})   );
      h_NOMCH_gen_mu_dEta01_2->Fill( abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(1).Eta() ) );
      h_NOMCH_gen_mu_dEta02_2->Fill( abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(2).Eta() ) );
      h_NOMCH_gen_mu_dEta12_2->Fill( abs( Muons.at(i).at(1).Eta() - Muons.at(i).at(2).Eta() ) );
      h_NOMCH_gen_mu_dEtamin_2->Fill(std::min( {abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(1).Eta() ) , abs( Muons.at(i).at(0).Eta() - Muons.at(i).at(2).Eta() ) , abs( Muons.at(i).at(1).Eta() - Muons.at(i).at(2).Eta() )   })   );
      h_NOMCH_gen_mu_pt_eta_2->Fill( Muons.at(i).at(2).Pt(), abs(Muons.at(i).at(2).Eta()  ) );
      h_NOMCH_gen_mu_n_eta_2->Fill( Tracks.at(i).size(), abs(Muons.at(i).at(2).Eta()  ) );
      h_NOMCH_gen_mu_pt_n_2->Fill( Muons.at(i).at(2).Pt(), Tracks.at(i).size() );
    }
    h_MCH_nmatches->Fill(matches);

    if(matches==3)
    {
    h_MCH_invmass->Fill(    ( get<0>(get<1>(flags.at(0)))   +  get<0>(get<1>(flags.at(1)) ) +  get<0>(get<1>(flags.at(2))  )     ).M()   );  
    }
  }

}

std::vector<std::tuple<bool,std::tuple<TLorentzVector,int,float>>> RunMatching(std::vector<TLorentzVector> Muons,std::vector<std::tuple<TLorentzVector,int,float>> Tracks)
{
  //define indexes vector
  int indexes[3] = {-1,-1,-1};
  //Get index of matched track to the three muons
  //Start from the third muon
  for(uint j=2; j<3; ++j)
  {
    int index=-1;float w = 0.01;
    for(uint k=0; k<Tracks.size(); ++k)
    {
     float dR = Muons.at(j).DeltaR( get<0>(Tracks.at(k)));
     if(dR<w && !FindValue(k,indexes) ){w=dR;index=k;}
    }
    indexes[j]=index; 
  }
  //Start from the second muon
  for(uint j=1; j<2; ++j)
  {
    int index=-1;float w = 0.01;
    for(uint k=0; k<Tracks.size(); ++k)
    {
     float dR = Muons.at(j).DeltaR(get<0>(Tracks.at(k)));
     if(dR<w && !FindValue(k,indexes) ){w=dR;index=k;}
    }
    indexes[j]=index; 
  }
  //Start from the first muon
  for(uint j=0; j<1; ++j)
  {
    int index=-1;float w = 0.01;
    for(uint k=0; k<Tracks.size(); ++k)
    {
     float dR = Muons.at(j).DeltaR( get<0>(Tracks.at(k))  );
     if(dR<w && !FindValue(k,indexes) ){w=dR;index=k;}
    }
    indexes[j]=index; 
  }

  //Make output
  std::vector<std::tuple<bool,std::tuple<TLorentzVector,int,float>>> output;
  for(uint j=0; j<Muons.size(); ++j)
  {
    int value=indexes[j];
    std::tuple<TLorentzVector,int,float> track;
    if( value != -1 )
    {
     track = Tracks.at(value);
     output.push_back( make_tuple(true,track) );
    } 
    else
    {
     output.push_back( make_tuple(false,track) );
    }
  }
  return output;
}

bool FindValue(int index, int values[3])
{
    bool found=false;
    if(values[0]==index) found=true;
    if(values[1]==index) found=true;
    if(values[2]==index) found=true;
    return found;
}