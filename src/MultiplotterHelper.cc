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
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TLine.h>
#include "/uscms/home/guerrero/nobackup/Phase2/tau3mufwk/scripts/CMS_lumi.c"
#include "MultiplotterHelper.h"

using namespace std;


void MakeRatePlot(TH1F* h_variable, float scale, string tagname)
{
//Create rate histogram
TH1F *h_rate = new TH1F(
    (std::string(h_variable->GetName())+std::string( Form("_%s",tagname.c_str() ))).c_str(),
    (std::string("rate")+std::string(h_variable->GetTitle())).c_str(),
    h_variable->GetNbinsX(),
    h_variable->GetBinLowEdge(1),
    h_variable->GetBinLowEdge(h_variable->GetNbinsX()+1)
);

//Rate (kHz) Calculating and scaling
float total_rate    = scale;
float total_entries = h_variable->Integral(-1, -1); 

for (int ibin = 1; ibin < h_variable->GetNbinsX()+1; ++ibin)
{
    float tmp_integral   = h_variable->Integral(ibin, -1);
    float rate_reduction = tmp_integral/total_entries;
    h_rate->SetBinContent(ibin, rate_reduction);
}
h_rate->Scale(total_rate);
}

TCanvas* CreateAcceptanceMapPlot(TH2F *tmp, string name)
{
   TCanvas *c1 = new TCanvas(Form("%s",name.c_str()),Form("%s",name.c_str()),10,10,1000,1000);
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0000);
   gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
   gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
   gStyle->SetLineWidth(3);
   c1->SetGrid();
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);
   tmp->Draw("COLZ TEXT");
   return c1;
}

vector< tuple<TGraphAsymmErrors*, string> > TGraphProducer(string filename, string historef,const vector<string>& histonames, const vector<string>& legends){

  //Open file and ref histogram
  TFile *f_data=new TFile(filename.c_str());
  TH1F *h_ref=(TH1F*)f_data->Get( Form("%s",historef.c_str()));
  //Read other histograms
  int nHistos= histonames.size();
  TH1F *myhisto[nHistos] = {NULL};
  for(int histogramIT=0; histogramIT<nHistos; histogramIT++){
      bool exists = f_data->GetListOfKeys()->Contains(Form("%s",histonames[histogramIT].c_str()));
      if(!exists) {cout<<"Histogram "<<histonames[histogramIT].c_str()<<"does not exist on"<<filename<<endl;}
      if(!exists) continue;  
      myhisto[histogramIT] =(TH1F*)f_data->Get(Form("%s",histonames[histogramIT].c_str())); 
  } 
  //Create a vector of tgraphs
  vector< tuple<TGraphAsymmErrors*, string> > output;
  string tmpstr;
  for(int histogramIT=0; histogramIT<nHistos; histogramIT++){ 
      tmpstr = legends[histogramIT];
      TGraphAsymmErrors* tmp = new TGraphAsymmErrors(myhisto[histogramIT],h_ref); 
      output.push_back( make_tuple(tmp,tmpstr)  ); 
  }   
  return output;
}

void MultiEfficiencyCurves( vector< tuple<TGraphAsymmErrors*, string> > curves, const vector<int>& colors, string xlabel, string ylabel, string samplename, string pileup, const vector<string>& cutlegends, string plotname){
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0000);
  gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(3);
  int iPeriod = 5;     
  int iPos = 11;
  writeExtraText = true;       
  extraText  = "Phase-2 Simulation"; 
  lumi_13TeV = Form("%s, %s PU",samplename.c_str(),pileup.c_str());
  

  //Get the curves and asssign them a plot style
  int nCurves= curves.size();
  TGraphAsymmErrors *myCurves[nCurves] = {NULL};
  for(int k=0; k<nCurves; k++){ 
    myCurves[k] = get<0>(curves[k]);
    myCurves[k]->SetMarkerStyle(kFullSquare);
    myCurves[k]->SetLineWidth(4);
    myCurves[k]->SetLineColor(colors.at(k));
    myCurves[k]->SetMarkerColor(colors.at(k));
    myCurves[k]->SetMarkerSize(2);
  } 

  //Create canvas
  TCanvas *c_d=new TCanvas("c_d", "c_d", 800, 800);
  c_d->SetFillStyle(4000); c_d->SetFrameFillColor(0);
  TPad *p_d=new TPad("p_d", "p_d", 0, 0, 1, 1);
  p_d->SetFillStyle(4000);p_d->SetFrameFillColor(0);
  p_d->SetTopMargin(0.10);p_d->SetBottomMargin(0.10);
  p_d->SetLeftMargin(0.16);p_d->SetRightMargin(0.05);           
  p_d->Draw();     
  p_d->cd();  
  //Draw main curve
  myCurves[0]->Draw();
  //Draw the other curves
  for(uint i=1;i<curves.size();i++){myCurves[i]->Draw("SAME lp");}  
  myCurves[0]->SetMaximum(1.5);  
  myCurves[0]->SetMinimum(0);   
  //myCurves[0]->GetXaxis()->SetRangeUser(40,400);
  myCurves[0]->GetXaxis()->SetLabelFont(42);
  myCurves[0]->GetYaxis()->SetLabelFont(42);
  myCurves[0]->GetYaxis()->SetTitleFont(42);
  myCurves[0]->GetXaxis()->SetTitleFont(42);     
  myCurves[0]->GetYaxis()->SetTitleOffset(1.5);
  myCurves[0]->GetYaxis()->SetTitleSize(0.045);
  myCurves[0]->GetYaxis()->SetLabelSize(0.040);
  myCurves[0]->GetXaxis()->SetTitleSize(0.045);
  myCurves[0]->GetXaxis()->SetLabelSize(0.040);
  myCurves[0]->SetTitle( Form(";%s; %s", xlabel.c_str(),ylabel.c_str() )   ); 

  TLegend *leg = new TLegend(0.49,0.65,0.83,0.85,NULL,"brNDC");
  leg->SetNColumns(1);leg->SetBorderSize(0);
  leg->SetTextSize(0.023);leg->SetTextFont(42);
  leg->SetLineColor(1);leg->SetLineWidth(10);
  leg->SetFillColor(0);leg->SetFillStyle(0);
  for(uint j=0;j<curves.size();j++){leg->AddEntry(myCurves[j],Form("%s",get<1>(curves[j]).c_str() ), "lp");}
  leg->Draw();

  TPaveText *t1 = new TPaveText(0.20,0.75,0.48,0.88,"brNDC");
  t1->SetTextAlign(22);t1->SetTextFont(42);
  t1->SetTextSize(0.03);t1->SetFillColor(kWhite);
  t1->SetFillStyle(1001);t1->SetLineColor(1);
  t1->SetLineWidth(10);t1->SetBorderSize(0);
  for (uint cutsIT=0; cutsIT<cutlegends.size(); cutsIT++){
     t1->AddText(Form("%s", cutlegends[cutsIT].c_str()));
  }  
  t1->Draw();

  CMS_lumi( p_d, iPeriod, iPos );  
  c_d->SaveAs(Form("h_%s.png",plotname.c_str() ));

  delete c_d;

}

void MultiCanvas2D(string filename, const vector<string>& histonames, string xlabel, string ylabel, string samplename, string pileup, const vector<string>& cutlegends)
{

  //Style
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0000);
  gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(2);  
  gStyle->SetPalette(55);
  int iPeriod = 5;     
  int iPos = 11;
  writeExtraText = true;       
  extraText  = "Phase-2 Simulation"; 
  lumi_13TeV = Form("%s, %s PU",samplename.c_str(),pileup.c_str());
  TFile *f_data=new TFile(filename.c_str());
  int nHistos= histonames.size();
  TH2F *myhisto[nHistos] = {NULL};
  TCanvas *mycanvas[nHistos] = {NULL};
  for(int histogramIT=0; histogramIT<nHistos; histogramIT++){
      bool exists = f_data->GetListOfKeys()->Contains(Form("%s",histonames[histogramIT].c_str()));
      if(!exists) {cout<<"Histogram "<<histonames[histogramIT].c_str()<<"does not exist on"<<filename<<endl;}
      if(!exists) continue;
      myhisto[histogramIT] =(TH2F*)f_data->Get(Form("%s",histonames[histogramIT].c_str()));
      mycanvas[histogramIT]=new TCanvas(Form("Canvas_%i",histogramIT),Form("Canvas_%i",histogramIT) , 2500, 2000);
      mycanvas[histogramIT]->SetFillStyle(4000); mycanvas[histogramIT]->SetFrameFillColor(0);
      TPad *p_d=new TPad("p_d", "p_d", 0, 0, 1, 1);
      p_d->SetFillStyle(4000);p_d->SetFrameFillColor(0);
      p_d->SetBottomMargin(0.10);p_d->SetLeftMargin(0.18);p_d->SetRightMargin(0.15);           
      p_d->Draw();     
      p_d->cd();   
      myhisto[histogramIT]->Draw("COLZ");
      myhisto[histogramIT]->GetXaxis()->SetLabelFont(42); myhisto[histogramIT]->GetYaxis()->SetTitleSize(0.06);
      myhisto[histogramIT]->GetYaxis()->SetLabelFont(42); myhisto[histogramIT]->GetYaxis()->SetLabelSize(0.05);
      myhisto[histogramIT]->GetYaxis()->SetTitleFont(42); myhisto[histogramIT]->GetXaxis()->SetTitleSize(0.05);
      myhisto[histogramIT]->GetXaxis()->SetTitleFont(42); myhisto[histogramIT]->GetXaxis()->SetLabelSize(0.05); 
      myhisto[histogramIT]->GetYaxis()->SetTitleOffset(1.5);
      myhisto[histogramIT]->SetTitle(Form(";%s;%s",xlabel.c_str(),ylabel.c_str()) ); 
      p_d->RedrawAxis("g");

      TPaveText *t1 = new TPaveText(0.45,0.77,0.80,0.85,"brNDC");
      t1->SetTextAlign(22);t1->SetTextFont(42);
      t1->SetTextSize(0.03);t1->SetFillColor(kWhite);
      t1->SetFillStyle(1001);t1->SetLineColor(1);
      t1->SetLineWidth(10);t1->SetBorderSize(0);
      for (uint cutsIT=0; cutsIT<cutlegends.size(); cutsIT++){
         if (cutlegends[cutsIT].substr(0,4) !="None")
         {
         t1->AddText(Form("%s", cutlegends[cutsIT].c_str()));
         } 
      }  
      t1->Draw();

      CMS_lumi( p_d, iPeriod, iPos );      
      mycanvas[histogramIT]->SaveAs(Form("%s.png",histonames[histogramIT].c_str() ) );
      delete mycanvas[histogramIT];  
  }
}

void MultiCanvas1D(string filename, const vector<string>& histonames, 
  const vector<string>& histolegends,
  const vector<int>& histocolors, string xlabel, string ylabel, float maxyaxis, float minyaxis,
  int logyscale, int norm, string samplename, string pileup, const vector<string>& cutlegends)
{

  gROOT->SetStyle("Plain");gStyle->SetOptStat(0000);
  gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(3);
  int iPeriod = 5;     
  int iPos = 11;
  writeExtraText = true;       
  extraText  = "Phase-2 Simulation"; 
  lumi_13TeV = Form("%s, %s PU",samplename.c_str(),pileup.c_str()); 

  TFile *f_data=new TFile(filename.c_str());
  int nHistos= histonames.size();
  TH1F *myhisto[nHistos] = {NULL};
  TCanvas *c_d=new TCanvas("c_d", "c_d", 800, 800);
  c_d->SetFillStyle(4000); c_d->SetFrameFillColor(0);
  TPad *p_d=new TPad("p_d", "p_d", 0, 0, 1, 1);
  p_d->SetFillStyle(4000);p_d->SetFrameFillColor(0);
  p_d->SetTopMargin(0.10);p_d->SetBottomMargin(0.10);
  p_d->SetLeftMargin(0.16);p_d->SetRightMargin(0.05);           
  p_d->Draw();     
  p_d->cd(); 
  //Log Y scale plot
  if(logyscale==1) gPad->SetLogy(); 
  //Check which one to plot first
  float first=0,max;
  int idx=0;
  for(int histogramIT=0; histogramIT<nHistos; histogramIT++){
      bool exists = f_data->GetListOfKeys()->Contains(Form("%s",histonames[histogramIT].c_str()));
      if(!exists) {cout<<"Histogram "<<histonames[histogramIT].c_str()<<"does not exist on"<<filename<<endl;}
      if(!exists) continue;  
      myhisto[histogramIT] =(TH1F*)f_data->Get(Form("%s",histonames[histogramIT].c_str()));       
      if(norm==1) myhisto[histogramIT]->Scale(1. / myhisto[histogramIT]->Integral() ); 
      max = myhisto[histogramIT]->GetBinContent(myhisto[histogramIT]->GetMaximumBin());
      if(max>first){first=max;idx=histogramIT;}
  } 
  
  //Draw histogram with highest bin first  
  myhisto[idx]->SetLineWidth(4);   
  myhisto[idx]->SetLineColor(histocolors[idx]); 
  myhisto[idx]->GetXaxis()->SetLabelFont(42); myhisto[idx]->GetYaxis()->SetTitleSize(0.045);
  myhisto[idx]->GetYaxis()->SetLabelFont(42); myhisto[idx]->GetYaxis()->SetLabelSize(0.045);
  myhisto[idx]->GetYaxis()->SetTitleFont(42); myhisto[idx]->GetXaxis()->SetTitleSize(0.045);
  myhisto[idx]->GetXaxis()->SetTitleFont(42); myhisto[idx]->GetXaxis()->SetLabelSize(0.045);      
  myhisto[idx]->GetYaxis()->SetTitleOffset(1.75);
  myhisto[idx]->SetTitle(Form(";%s;%s",xlabel.c_str(),ylabel.c_str()) ); 
  if(logyscale==1) 
  {
    myhisto[idx]->Draw("HISTO");myhisto[idx]->SetMaximum(maxyaxis);myhisto[idx]->SetMinimum(minyaxis);
  }
  else
  {
    maxyaxis=1.5*first;
    myhisto[idx]->Draw("HISTO");myhisto[idx]->SetMaximum(maxyaxis);
  }
  //Draw the rest
  for(int histogramIT=0; histogramIT<nHistos; histogramIT++){ 
      bool exist = f_data->GetListOfKeys()->Contains(Form("%s",histonames[histogramIT].c_str()));
      if(!exist) {cout<<"Histogram "<<histonames[histogramIT].c_str()<<"does not exist on"<<filename<<endl;}
      if(!exist) continue;   
      if(histogramIT==idx) continue;  
      myhisto[histogramIT]->SetLineWidth(4);  
      myhisto[histogramIT]->SetLineColor(histocolors[histogramIT]); 
      if(norm==1){myhisto[histogramIT]->DrawNormalized("Histo SAME");}
      else{myhisto[histogramIT]->Draw("Histo SAME");}  
  }   

  //Draw  histolegends & cutlegends
  if(cutlegends.size() != 0)
  {
      if (histolegends[0].substr(0,4) !="None")
      {
      TLegend *leg_d = new TLegend(0.49,0.70,0.85,0.88,NULL,"brNDC");
      leg_d->SetNColumns(1);leg_d->SetBorderSize(0);
      leg_d->SetTextSize(0.023);leg_d->SetTextFont(42);
      leg_d->SetLineColor(1);leg_d->SetLineWidth(10);
      leg_d->SetFillColor(0);leg_d->SetFillStyle(0);
      for(int histogramIT=0; histogramIT<nHistos; histogramIT++){leg_d->AddEntry(myhisto[histogramIT], Form("%s", histolegends[histogramIT].c_str()) , "lp");}  
      leg_d->Draw();
      } 
      TPaveText *t1 = new TPaveText(0.20,0.75,0.48,0.88,"brNDC");
      t1->SetTextAlign(22);t1->SetTextFont(42);
      t1->SetTextSize(0.03);t1->SetFillColor(kWhite);
      t1->SetFillStyle(4000);t1->SetLineColor(1);
      t1->SetLineWidth(10);t1->SetBorderSize(0);
      for (uint cutsIT=0; cutsIT<cutlegends.size(); cutsIT++){
         t1->AddText(Form("%s", cutlegends[cutsIT].c_str()));
      }  
      t1->Draw();
  }
  else
  {
      if (histolegends[0].substr(0,4) !="None")
      {
      TLegend *leg_d = new TLegend(0.19,0.70,0.85,0.88,NULL,"brNDC");
      leg_d->SetNColumns(1);leg_d->SetBorderSize(0);
      leg_d->SetTextSize(0.03);leg_d->SetTextFont(42);
      leg_d->SetLineColor(1);leg_d->SetLineWidth(10);
      leg_d->SetFillColor(0);leg_d->SetFillStyle(0);
      for(int histogramIT=0; histogramIT<nHistos; histogramIT++){leg_d->AddEntry(myhisto[histogramIT], Form("%s", histolegends[histogramIT].c_str()) , "lp");}  
      leg_d->Draw();
      }
  } 

  CMS_lumi( p_d, iPeriod, iPos );
  p_d->RedrawAxis("g");  
  if(norm==1){
  if(nHistos>1){c_d->SaveAs(Form("%s_more_normalized.png",histonames[0].c_str() ) ); }
  else{c_d->SaveAs(Form("%s_normalized.png",histonames[0].c_str() ) ); }
  }
  else{
  if(nHistos>1){c_d->SaveAs(Form("%s_more.png",histonames[0].c_str() ) ); }
  else{c_d->SaveAs(Form("%s.png",histonames[0].c_str() ) ); }    
  }

  delete c_d;
}