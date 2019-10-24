#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <algorithm>
#include <cassert>
#include "CfgParser.h"
#include <tuple>
#include <iterator>
#include <vector>
#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TColor.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include "MultiplotterHelper.h"

using namespace std;

void PlotsEff(string filename, CfgParser configEff,string samplename, string pileup)
{
    //Parse over the list of sets
    std::vector<std::string> setEff;
    for (auto set : configEff.readStringListOpt("Booking::List"))
             setEff.push_back(set); 
    int nSet= setEff.size();
    std::vector<std::string> cutlegends;
    for (auto cutlegend : configEff.readStringListOpt("CutLegends::List"))
             cutlegends.push_back(cutlegend);
    //Loop over the defined sets
    for(int setIT=0; setIT<nSet; setIT++){  
    std::vector<std::string>histonames,legends;
    std::vector<int> colors;
    string xlabel,ylabel,historefname,plotname;
    vector< tuple<TGraphAsymmErrors*, string> > graphs;
        //Parse over the histo names/histo legend/histo color
        for (auto name : configEff.readStringListOpt(Form("%s::Historef", setEff[setIT].c_str())))
             historefname = name;  
        for (auto name : configEff.readStringListOpt(Form("%s::Histos", setEff[setIT].c_str())))
             histonames.push_back(name);      
        for (auto legend : configEff.readStringListOpt(Form("%s::Legends", setEff[setIT].c_str())))
             legends.push_back(legend);
        for (auto color : configEff.readStringListOpt(Form("%s::Colors", setEff[setIT].c_str() )))
             colors.push_back(stoi(color));
        for (auto x : configEff.readStringListOpt(Form("%s::Xlabel", setEff[setIT].c_str() )) )
             xlabel = x;
        for (auto y : configEff.readStringListOpt(Form("%s::Ylabel", setEff[setIT].c_str() )) )
             ylabel = y;
        for (auto p : configEff.readStringListOpt(Form("%s::Plotname", setEff[setIT].c_str() )) )
             plotname = p;
        graphs = TGraphProducer(filename,historefname,histonames,legends);
        MultiEfficiencyCurves(graphs,colors,xlabel,ylabel,samplename,pileup,cutlegends,plotname);
    }
}

void Plots1D(string filename, CfgParser config1D, string samplename, string pileup)
{
    //Parse over the list of sets
    std::vector<std::string> set1D;
    for (auto set : config1D.readStringListOpt("Booking::List"))
             set1D.push_back(set); 
    int nSet= set1D.size();
    std::vector<std::string> cutlegends;
    for (auto cutlegend : config1D.readStringListOpt("CutLegends::List"))
             cutlegends.push_back(cutlegend);
    //Loop over the defined sets
    for(int setIT=0; setIT<nSet; setIT++){  
    std::vector<std::string>histonames,histolegends;
    std::vector<int> histocolors;
    int norm,logyscale;
    float minyaxis,maxyaxis;    
    string xlabel,ylabel;
        //Parse over the histo names/histo legend/histo color
        for (auto name : config1D.readStringListOpt(Form("%s::Histos", set1D[setIT].c_str())))
             histonames.push_back(name);      
        for (auto legend : config1D.readStringListOpt(Form("%s::Legend", set1D[setIT].c_str())))
             histolegends.push_back(legend);
        for (auto colors : config1D.readStringListOpt(Form("%s::Colors", set1D[setIT].c_str() )))
             histocolors.push_back(stoi(colors));
        for (auto x : config1D.readStringListOpt(Form("%s::Xlabel", set1D[setIT].c_str() )) )
             xlabel = x;
        for (auto y : config1D.readStringListOpt(Form("%s::Ylabel", set1D[setIT].c_str() )) )
             ylabel = y;  
         for (auto n : config1D.readStringListOpt(Form("%s::MaxYaxis", set1D[setIT].c_str() )) )
             maxyaxis = stof(n);
        for (auto n : config1D.readStringListOpt(Form("%s::MinYaxis", set1D[setIT].c_str() )) )
             minyaxis = stof(n);
        for (auto n : config1D.readStringListOpt(Form("%s::LogYscale", set1D[setIT].c_str() )) )
             logyscale = stoi(n); 
        for (auto n : config1D.readStringListOpt(Form("%s::Norm", set1D[setIT].c_str() )) )
             norm = stoi(n); 
        MultiCanvas1D(filename,histonames,histolegends,histocolors,xlabel,ylabel,maxyaxis,minyaxis,logyscale,norm,samplename,pileup,cutlegends);
    }
}

void Plots2D(string filename, CfgParser config2D, string samplename, string pileup)
{
    //Parse over the list of sets
    std::vector<std::string> set2D;
    for (auto set : config2D.readStringListOpt("Booking::List"))
             set2D.push_back(set); 
    int nSet= set2D.size();
    std::vector<std::string> cutlegends;
    for (auto cutlegend : config2D.readStringListOpt("CutLegends::List"))
             cutlegends.push_back(cutlegend);
    //Loop over the defined sets
    for(int setIT=0; setIT<nSet; setIT++){  
    std::vector<std::string>histonames;
    string xlabel,ylabel;
        //Parse over the histo names/histo legend/histo color
        for (auto name : config2D.readStringListOpt(Form("%s::Histos", set2D[setIT].c_str())))
             histonames.push_back(name);      
        for (auto x : config2D.readStringListOpt(Form("%s::Xlabel", set2D[setIT].c_str() )) )
             xlabel = x;
        for (auto y : config2D.readStringListOpt(Form("%s::Ylabel", set2D[setIT].c_str() )) )
             ylabel = y;        
        MultiCanvas2D(filename,histonames,xlabel,ylabel,samplename,pileup,cutlegends); 
    }
}

const bool filecheck(string filename){
    if ( access(filename.c_str(), F_OK ) == -1 ){
        cout<<"[WARNING] "<<filename<<" file does not exist or it was not specified"<<endl;
        return false;        
    }
    else{
        return true;
    }
}

int main(int argc, char** argv)
{
    cout << "[INFO] ... starting program to plot some physics" << endl;
    
    po::options_description desc("Plotting options");
    desc.add_options()
        ("help", "produce help message")
        ("inputFile" , po::value<string>()->required(), "input root file name")       
        ("samplename" , po::value<string>()->required(), "Sample name") 
        ("pileup" , po::value<string>()->required(), "Pile up scenario")    
        //Optional    
        ("configEff" , po::value<string>()->default_value("cfg"), "name of configuration file with efficiency-graph information (e.g. config/*.cfg)")
        ("configTH1" , po::value<string>()->default_value("cfg"), "name of configuration file with TH1-histogram information (e.g. config/*.cfg)")
        ("configTH2" , po::value<string>()->default_value("cfg"), "name of configuration file with TH2-histogram information (e.g. config/*.cfg)")
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
    
    cout << "[INFO] ... searching for root file with histograms " << opts["inputFile"].as<string>().c_str() << endl;
    const bool infoFile   = filecheck(opts["inputFile"].as<string>().c_str());
    string samplename = opts["samplename"].as<string>().c_str();
    string pileup     = opts["pileup"].as<string>().c_str();
    if(!infoFile){cerr << "** [ERROR] Aborting ..." << endl;return 1;}  

    //1D Plots 
    const bool info1D = filecheck(opts["configTH1"].as<string>().c_str());     
    if (info1D)
    {
        //Make 1D Plots in a canvas
        CfgParser config1D;
        config1D.init(opts["configTH1"].as<string>()) ;
        Plots1D(opts["inputFile"].as<string>().c_str(), config1D, samplename, pileup);
        //Save it in directory defined in config file
        string outputFolder1D = config1D.readStringOpt("General::outputFolder");        
        system ("mkdir myplots");
        system (Form("mkdir myplots/%s_%s_%s", samplename.c_str(), pileup.c_str(), outputFolder1D.c_str() ));
        system (Form("mv *.png *.pdf myplots/%s_%s_%s",samplename.c_str(), pileup.c_str(), outputFolder1D.c_str() ));       
    }
    else{cout<<"[WARNING] No information for 1D histogram plots"<<endl;}

    //2D Plots
    const bool info2D = filecheck(opts["configTH2"].as<string>().c_str());
    if (info2D)
    {
        //Make 2D Plots in a canvas
        CfgParser config2D;
        config2D.init(opts["configTH2"].as<string>()) ;
        Plots2D(opts["inputFile"].as<string>().c_str(), config2D,samplename, pileup);        
        //Save it in directory defined in config file
        string outputFolder2D = config2D.readStringOpt("General::outputFolder");
        system ("mkdir myplots");
        system (Form("mkdir myplots/%s_%s_%s", samplename.c_str(), pileup.c_str(), outputFolder2D.c_str() ));
        system (Form("mv *.png *.pdf myplots/%s_%s_%s",samplename.c_str(), pileup.c_str(), outputFolder2D.c_str() ));    
    }
    else{cout<<"[WARNING] No information for 2D histogram plots"<<endl;}


    //Efficiencys Plots 
    const bool infoEff = filecheck(opts["configEff"].as<string>().c_str());     
    if (infoEff)
    {
        //Make 1D Plots in a canvas
        CfgParser configEff;
        configEff.init(opts["configEff"].as<string>()) ;
        PlotsEff(opts["inputFile"].as<string>().c_str(), configEff,samplename, pileup);
        //Save it in directory defined in config file
        string outputFolderEff = configEff.readStringOpt("General::outputFolder");        
        system ("mkdir myplots");     
        system (Form("mkdir myplots/%s_%s_%s", samplename.c_str(), pileup.c_str(), outputFolderEff.c_str() ));
        system (Form("mv *.png *.pdf myplots/%s_%s_%s",samplename.c_str(), pileup.c_str(), outputFolderEff.c_str() ));
    }
    else{cout<<"[WARNING] No information for efficiency plots"<<endl;}

}
