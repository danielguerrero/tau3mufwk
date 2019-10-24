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

void MuonTriggerEfficiencies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons, 
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons, 
	string outputname);

void TriggerEfficiencyStudies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons);

void PureMuonTriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks, float mwindow, float drwindow, float dzwindow, string tagname);
void MixedMuonTriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Track3, float mwindow, float drwindow, float dzwindow, string tagname);
void MixedME0TriggerCalculator(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Tracks12,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> Track3, float dzwindow, float drwindow, float bendcut, string tagname);


void TriggerOREfficiencyCalculator(std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons, int trig1, int trig2, int trig3);