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

void MuonEfficiencies(std::vector<std::vector<TLorentzVector>> GenMuons,
	std::vector<std::vector<TLorentzVector>> EMTFMuons,
	std::vector<std::vector<TLorentzVector>> L1TTMuons,
	std::vector<std::vector<TLorentzVector>> L1TkMuMuons,
	std::vector<std::vector<TLorentzVector>> L1TkMuStubMuons,string outputname);

void EfficienciesStudies(std::vector<std::vector<TLorentzVector>> genmuons, std::vector<std::vector<TLorentzVector>> tracks, string tagname);