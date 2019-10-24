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

void TriggerRateStudies(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons);

void MuonTriggerRates(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuons,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuons, 
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuons, 	
	string outputname);

void TriggerORRateCalculator(std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TTMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> L1TkMuStubMuons,
  std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> ME0StubMuons, int trig1, int trig2, int trig3);