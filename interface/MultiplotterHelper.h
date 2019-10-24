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

using namespace std;

void MakeRatePlot(TH1F* h_variable, float scale, string tagname);

TCanvas* CreateAcceptanceMapPlot(TH2F *tmp, string name);

vector< tuple<TGraphAsymmErrors*, string> > TGraphProducer(string filename, string historef,const vector<string>& histonames, const vector<string>& legends);

void MultiEfficiencyCurves( vector< tuple<TGraphAsymmErrors*, string> > curves, const vector<int>& colors, string xlabel, string ylabel,string sample, string pileup, const vector<string>& cutlegends, string plotname);

void MultiCanvas2D(string filename, const vector<string>& histonames, string xlabel, string ylabel, string sample, string pileup, const vector<string>& cutlegends);

void MultiCanvas1D(string filename, const vector<string>& histonames, 
  const vector<string>& histolegends,
  const vector<int>& histocolors, string xlabel, string ylabel, float maxyaxis, float minyaxis, int logyscale, int norm, string sample, string pileup,const vector<string>& cutlegends);