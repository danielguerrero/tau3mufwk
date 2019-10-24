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
#include "CfgParser.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>  
#include "SelectionHelper.h"
#include "GenProperties.h"
#include "MuonTrackProperties.h"
#include "MuonEfficiencies.h"
#include "MuonHits.h"
#include "MuonTriggerEfficiencies.h"
#include "MuonTriggerRates.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

const bool filecheck(string filename){
	if ( access(filename.c_str(), F_OK ) == -1 ){
		cout<<"[ERROR] "<<filename<<" file does not exist or it was not specified"<<endl;
		return false;        
	}
	else{
		return true;
	}
}


int main (int argc, char** argv) {
	cout << "[INFO] ... starting program" << endl;
	////////////////////////////////////////////////////////////////////////
	// Declare command line options
	////////////////////////////////////////////////////////////////////////
	po::options_description desc("Skim options");
	desc.add_options()
		("help", "produce help message")
		// required
		("input" , po::value<string>()->required(), "input file list or root file")
		("tagname", po::value<string>()->required(), "output tagname")
		("is-signal",   po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a tau->3mu signal sample (default is false)")
		("is-minbias",  po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a min bias sample (default is false)")
		("config" , po::value<string>()->default_value("cfg"), "name of configuration file with kinematic information (e.g. config/*.cfg)")
	;

	po::variables_map opts;
	try {
		po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), opts);
		if (opts.count("help")) {
			cout << desc << "\n";
			return 1;
		}
		po::notify(opts);
	}    
	catch (po::error& e) {
		cerr << "** [ERROR] " << e.what() << endl;
		return 1;
	}
	const bool is_signal  = opts["is-signal"].as<bool>();
	const bool is_minbias = opts["is-minbias"].as<bool>();
	////////////////////////////////////////////////////////////////////////
	// Prepare event loop
	////////////////////////////////////////////////////////////////////////
	cout << "[INFO] ... opening file list : " << opts["input"].as<string>().c_str() << endl;
	if ( access( opts["input"].as<string>().c_str(), F_OK ) == -1 ){
		cerr << "** [ERROR] The input file list does not exist, aborting" << endl;
		return 1;        
	}
	TFile *file = new TFile(opts["input"].as<string>().c_str()); 
	TTree *tree = (TTree*)file->Get("Ntuplizer/MuonTrackTree");
	////////////////////////////////////////////////////////////////////////
	// Make studies, output histograms using config file parameters
	////////////////////////////////////////////////////////////////////////
	string outputFileName = opts["tagname"].as<string>();
	system("mkdir histograms");
	//cout << "[INFO] ... saving output to file : " << outputFileName << endl; 
	const bool info = filecheck(opts["config"].as<string>().c_str());     
	if (info)
	{
		//Read Event Selection Config file
		CfgParser config; config.init(opts["config"].as<string>());

		//Select Events using Config file and Return Objects (To be replaced by a map, tuple for now)
		//------------------- (0) GenTaus, (1) GenTaus, (2) EMTF++tracks, (3) L1TTtracks, (4) L1TkMutracks, (5) L1TkMuStubtracks, (6) Hits 
		// Signal then we do studies in the positive endcap/if min bias sample, then looks for tracks in either endcap 
		std::vector<std::tuple<std::vector<TLorentzVector>,std::vector<TLorentzVector>,
		std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,
		std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::tuple<TLorentzVector,int,float>>,
		std::vector<std::tuple<TLorentzVector,int,float>>,std::vector<std::array<float,11>> >> ParticlesTracksHits = SelectParticlesTracksHits(tree,config,is_signal);
		std::vector<std::vector<TLorentzVector>> GenTaus, GenMuons;
		std::vector<std::vector< std::tuple<TLorentzVector, int,float> >> EMTFs,L1TTs, L1TkMus, L1TkMuStubs, ME0Stubs;
		std::vector<std::vector<std::array<float,11>>> Hits;
		GenTaus     = GetObjectGens(ParticlesTracksHits,0); //If it's not signal is empty 
		GenMuons    = GetObjectGens(ParticlesTracksHits,1); //If it's not signal is empty 
		EMTFs       = GetObjectTracks(ParticlesTracksHits,0); 
		L1TTs       = GetObjectTracks(ParticlesTracksHits,1); 
		L1TkMus     = GetObjectTracks(ParticlesTracksHits,2); 
		L1TkMuStubs = GetObjectTracks(ParticlesTracksHits,3); 
		ME0Stubs    = GetObjectTracks(ParticlesTracksHits,4);
		Hits        = GetObjectHits(ParticlesTracksHits); 
		//Perform Data Analysis using objects and modules for generator level
		//------------------- GenMuons, GenTaus, EMTF++tracks, L1TTtracks, L1TkMutracks, L1TkMuStubtracks 
		if(is_signal)
		{
			GenProperties(GenTaus,GenMuons,outputFileName);
			MuonTriggerEfficiencies(GenMuons,EMTFs,L1TTs,L1TkMus,L1TkMuStubs,ME0Stubs,outputFileName);
		} 
		//Perform Data Analysis using objects and modules for l1 trigger level
		//------------------- EMTF++tracks, L1TTtracks, L1TkMutracks, L1TkMuStubtracks 
		MuonTrackProperties(EMTFs,L1TTs,L1TkMus,L1TkMuStubs,ME0Stubs,outputFileName);
		if (is_minbias)
		{ 
		MuonTriggerRates(EMTFs,L1TTs,L1TkMus,L1TkMuStubs,ME0Stubs,outputFileName);
		} 
			

	}
	else{cout<<"[ERROR] No config file for kinematic cuts"<<endl;}

}
