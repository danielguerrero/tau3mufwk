#include <iostream>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>
#include "scripts/CMS_lumi.c"

TGraph* ROCcurve(TH1F *h_SIG,TH1F *h_BKG){
  //Get info for tgraphs for ROC curve using BDT1
  const Int_t  n = h_SIG->GetNbinsX();
  Double_t SIG[n],BKG[n];

  for (Int_t i=1;i<=n;i++) {
    SIG[i] =h_SIG->GetBinContent(i);
    BKG[i] =1-h_BKG->GetBinContent(i);
    //if(h_SIG->GetBinContent(i)==h_SIG->GetBinContent(i+1)) break;
  } 
  TGraph *curve = new TGraph(n,SIG,BKG);
 
  return curve;
} 

void ROC1(int pu, int ext)
{
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0000);
  gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(3);  
  gStyle->SetErrorX(0);
  int iPeriod = 5;     
  int iPos = 11;
  writeExtraText = true;       
  extraText  = "Phase-2 Simulation"; 
  lumi_13TeV = ""; 

  //Input Files
  char* signal; 
  if(ext==0)
    signal = Form("histograms/Tau3muPU%i_fiducial_muontriggers.root",pu);
  else
    signal = Form("histograms/Tau3muPU%i_fiducialextended_muontriggers.root",pu);
  TFile *f_data1=new TFile(signal);
  TFile *f_data2=new TFile("histograms/MinBiasPU200_muontriggerrates.root");
  //Get HISTOGRAMS
  TH1F *h_SIG_TRG_1 =(TH1F*)f_data1->Get("h_TRG1_variable_roc");
  TH1F *h_SIG_TRG_2 =(TH1F*)f_data1->Get("h_TRG2_variable_roc");
  TH1F *h_SIG_TRG_3 =(TH1F*)f_data1->Get("h_TRG3_variable_roc");

  TH1F *h_BKG_TRG_1 =(TH1F*)f_data2->Get("h_TRG_1_variable_roc");
  TH1F *h_BKG_TRG_2 =(TH1F*)f_data2->Get("h_TRG_2_variable_roc");
  TH1F *h_BKG_TRG_3 =(TH1F*)f_data2->Get("h_TRG_3_variable_roc");

  //ROC Curves
  TGraph *Graph1 = ROCcurve(h_SIG_TRG_1,h_BKG_TRG_1);
  TGraph *Graph2 = ROCcurve(h_SIG_TRG_2,h_BKG_TRG_2);
  TGraph *Graph3 = ROCcurve(h_SIG_TRG_3,h_BKG_TRG_3);

  //Colors
  Graph1->SetLineColor(1);
  Graph1->SetLineWidth(4);
  Graph1->SetLineStyle(1);
  Graph1->SetMarkerColor(1);
  Graph1->SetMarkerStyle(21);
  Graph2->SetLineColor(2);
  Graph2->SetLineWidth(4);
  Graph2->SetLineStyle(1);
  Graph2->SetMarkerColor(2);
  Graph2->SetMarkerStyle(21);
  Graph3->SetLineColor(4);
  Graph3->SetLineWidth(4);
  Graph3->SetLineStyle(1);
  Graph3->SetMarkerColor(4);
  Graph3->SetMarkerStyle(21);

  TCanvas *c_1=new TCanvas("c_1", "c_1", 1600, 1400);
  c_1->SetFillStyle(0000);
  c_1->SetFrameFillColor(0);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 1);
  p_2->SetFillStyle(4000);p_2->SetFrameFillColor(0);
  p_2->SetTopMargin(0.10);p_2->SetBottomMargin(0.10);
  p_2->SetLeftMargin(0.16);p_2->SetRightMargin(0.05);  
  p_2->Draw();
  p_2->cd();
  //gPad->SetLogx(); 
  Graph1->SetTitle("");
  Graph1->GetXaxis()->SetTitle("#tau#rightarrow3#mu trigger efficiency");
  Graph1->GetYaxis()->SetTitle("MinBias background rejection");
  Graph1->GetXaxis()->SetLabelFont(42);
  Graph1->GetYaxis()->SetLabelFont(42);
  Graph1->GetYaxis()->SetTitleFont(42);
  Graph1->GetXaxis()->SetTitleFont(42);     
  Graph1->GetYaxis()->SetTitleOffset(1.5);
  Graph1->GetYaxis()->SetTitleSize(0.055);
  Graph1->GetYaxis()->SetLabelSize(0.040);
  Graph1->GetXaxis()->SetTitleSize(0.045);
  Graph1->GetXaxis()->SetLabelSize(0.040);
  Graph1->Draw("ALP");
  Graph1->GetXaxis()->SetRangeUser(0,0.4);
  Graph1->GetYaxis()->SetRangeUser(0.99967,1);
  Graph2->Draw("SAMELP");
  Graph3->Draw("SAMELP");

  TLine *line1=new TLine(0, 1-(10/(2760*11.246)), 0.4, 1-(10/(2760*11.246)));
  line1->SetLineWidth(7); line1->SetLineStyle(1); line1->SetLineColor(kOrange);
  line1->Draw("SAME");
  TLine *line2=new TLine(0, 1-(1/(2760*11.246)),  0.4, 1-(1/(2760*11.246)));
  line2->SetLineWidth(7); line2->SetLineStyle(2); line2->SetLineColor(kOrange);
  line2->Draw("SAME");
  c_1->Update();  
  TLegend *leg_1 = new TLegend(0.50,0.70,0.85,0.85,NULL,"brNDC");
  leg_1->SetNColumns(1);leg_1->SetBorderSize(0);
  leg_1->SetTextSize(0.03);leg_1->SetTextFont(42);
  leg_1->SetLineColor(0);leg_1->SetLineWidth(10);
  leg_1->SetFillColor(kWhite);leg_1->SetFillStyle(1001);
  leg_1->AddEntry(Graph1,"L1T_3Tk", "lp");  
  leg_1->AddEntry(Graph2,"L1T_3TkMuStub", "lp");  
  leg_1->AddEntry(Graph3,"L1T_3TkMu", "lp"); 
  leg_1->AddEntry(line2,"Bkg reduction at 1 kHz", "l");
  leg_1->AddEntry(line1,"Bkg reduction at 10 kHz", "l");
  leg_1->Draw();
  lumi_13TeV = ""; 
  CMS_lumi( p_2, iPeriod, iPos );  
  system ("mkdir myplots"); 
  system ("mkdir myplots/ROCcurves"); 
  char* output;
  if (ext==0)
     output = Form("myplots/ROCcurves/Tau3MuTriggersROC_same_%i.png",pu) ;
  else
     output = Form("myplots/ROCcurves/Tau3MuTriggersROC_same_%i_extended.png",pu) ;
  c_1->SaveAs(output);
}

