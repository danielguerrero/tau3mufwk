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
#include <vector>
#include <string.h>
#include <stdio.h>
#include "CfgParser.h" 
 
using namespace std;

//Tools to selection objects (gens,tracks,ME0stubs,hits)
std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>, std::vector<std::array<float,11>> >> SelectParticlesTracksHits(TTree *tree, CfgParser config, const bool is_signal);
std::vector<std::vector< std::tuple<TLorentzVector,int,float> >> GetObjectTracks(std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>, std::vector<std::array<float,11>> >> objects, int objectid);
std::vector<std::vector<TLorentzVector>>                         GetObjectGens(  std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>, std::vector<std::array<float,11>> >> objects, int objectid);
std::vector<std::vector<std::array<float,11>>>                   GetObjectHits(  std::vector< std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>, std::vector<std::array<float,11>> >> objects);
std::vector<std::tuple<TLorentzVector,int,float>> SelectTracks(vector<float>* pt, vector<float>* eta, vector<float>* phi, vector<float>* chi2, vector<int>* nstubs, vector<int>* charge, vector<float>* z,  bool isEMTF, bool isTT, bool overlapremoval, int endcap, CfgParser config);
std::vector<std::tuple<TLorentzVector,int,float>> SelectStubs(vector<float>* eta, vector<float>* phi, vector<float>* theta, vector<float>* bend,vector<int>* quality, vector<int>* type, int endcap);
std::vector<std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>>> SelectParticles(std::vector<TLorentzVector> taus, std::vector<TLorentzVector> muons, CfgParser config);
std::vector<std::array<float,11>> SelectHits(vector<float>* eta, vector<float>* phi, vector<float>* r_xy, vector<float>* z, vector<int>* type, vector<int>* station, vector<int>* ring, vector<int>* sector, vector<int>* neighbor, int endcap, CfgParser config);

//Other tools
std::vector<std::tuple<float,float,float>> RemoveOvelapppingTracks(std::vector<std::tuple<float,float,float>> info);
std::vector<TLorentzVector> BuildVectorOfTLorentzVector(vector<float>* P1, vector<float>* P2, vector<float>* P3,vector<float>* P4);
string createxlabel(string histoname);
std::vector<std::string> split(std::string str,std::string sep);