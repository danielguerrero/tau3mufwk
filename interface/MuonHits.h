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
void MuonHits(std::vector<std::vector<std::array<float,11>>> Hits, string outputname);

void HitPropertiesStudies(std::vector<std::vector<std::array<float,11>>> Hits);

void GetHitOccupancies(std::vector<std::vector<std::array<float,11>>> Hits, int type, int station, string tagname);

