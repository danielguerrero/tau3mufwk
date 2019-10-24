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
void MuonTrackProperties(std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> EMTFMuon,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TTMuon,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuMuon,
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> L1TkMuStubMuon, 
	std::vector<std::vector<std::tuple<TLorentzVector,int,float>>> ME0StubMuon, 
	string outputname);