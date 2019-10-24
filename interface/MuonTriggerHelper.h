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
#include <string>
#include "CfgParser.h"

using namespace std;

std::vector<std::tuple<TLorentzVector,int,float>> FindMixedTrackCandidates(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12,std::vector<std::tuple<TLorentzVector,int,float>> Tracks3);
std::vector<std::tuple<TLorentzVector,int,float>> MixedTriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12,std::vector<std::tuple<TLorentzVector,int,float>> Tracks3, float mwindow, float drwindow, float dzwindow);
std::vector<std::tuple<TLorentzVector,int,float>> MixedME0TriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks12,std::vector<std::tuple<TLorentzVector,int,float>> Tracks3, float dzwindow, float drwindow, float bendcut);

std::vector<std::tuple<TLorentzVector,int,float>> FindPureTrackCandidates(std::vector<std::tuple<TLorentzVector,int,float>> Tracks);
std::vector<std::tuple<TLorentzVector,int,float>> PureTriggerSelector(std::vector<std::tuple<TLorentzVector,int,float>> Tracks,float mwindow, float drwindow, float dzwindow);

void MixedTripleMuonTriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12, std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float mwindow, float drwindow, float dzwindow, string tagname);
void PureTripleMuonTriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>>  Tracks, float mwindow, float drwindow, float dzwindow, string tagname);
void MixedTripleME0TriggerCalculator(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12, std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks3, float dzwindow, float drwindow, float bendcut, string tagname);

//Matching tools
bool FindValue(int index, int values[3]);
std::vector<std::tuple<bool,std::tuple<TLorentzVector,int,float>>> RunMatching(std::vector<TLorentzVector> Muons,std::vector<std::tuple<TLorentzVector,int,float>> Tracks);
void MatchingCalculator(std::vector<std::vector<TLorentzVector>> Muons,std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks, string tagname);