#include <iostream>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TF1.h>
#include <TLine.h>
#include "/uscms/home/guerrero/nobackup/Phase2/tau3mufwk/scripts/CMS_lumi.c"

TGraph* ROCcurve(TH1F *h_SIG,TH1F *h_BKG){
  //Get info for tgraphs for ROC curve using BDT1
  const Int_t  n = h_SIG->GetNbinsX();
  Double_t SIG[n],BKG[n];
  for (Int_t i=1;i<n+1;i++) { 
    SIG[i] =h_SIG->GetBinContent(i);
    BKG[i] =1-h_BKG->GetBinContent(i);
    //if(h_SIG->GetBinContent(i)==h_SIG->GetBinContent(i+1)) break;
  } 
  TGraph *curve = new TGraph(n,SIG,BKG);
 
  return curve;
} 

void ROC2(int pu, int ext)
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
  lumi_13TeV = Form("%i PU",pu); 

  //Input Files
  char* signal; 
  if(ext==0)
    signal = Form("histograms/Tau3muPU%i_fiducial_muontriggers.root",pu);
  else
    signal = Form("histograms/Tau3muPU%i_fiducialextended_muontriggers.root",pu);
  TFile *f_data1=new TFile(signal);
  TFile *f_data2=new TFile("histograms/MinBiasPU200_muontriggerrates.root");
  //Get HISTOGRAMS
  TH1F *h_SIG_TRG_4  =(TH1F*)f_data1->Get("h_TRG4_variable_roc");
  TH1F *h_SIG_TRG_7  =(TH1F*)f_data1->Get("h_TRG5_variable_roc");
  TH1F *h_SIG_TRG_10 =(TH1F*)f_data1->Get("h_TRG6_variable_roc");
  TH1F *h_SIG_TRG_13 =(TH1F*)f_data1->Get("h_TRG7_variable_roc");
  TH1F *h_BKG_TRG_4 =(TH1F*)f_data2->Get("h_TRG_4_variable_roc");
  TH1F *h_BKG_TRG_7 =(TH1F*)f_data2->Get("h_TRG_5_variable_roc");
  TH1F *h_BKG_TRG_10 =(TH1F*)f_data2->Get("h_TRG_6_variable_roc");
  TH1F *h_BKG_TRG_13 =(TH1F*)f_data2->Get("h_TRG_7_variable_roc");

  //ROC Curves
  TGraph *Graph4 = ROCcurve(h_SIG_TRG_4,h_BKG_TRG_4);
  TGraph *Graph7 = ROCcurve(h_SIG_TRG_7,h_BKG_TRG_7);
  TGraph *Graph10 = ROCcurve(h_SIG_TRG_10,h_BKG_TRG_10); 
  TGraph *Graph13 = ROCcurve(h_SIG_TRG_13,h_BKG_TRG_13);

  Graph4->SetLineColor(6);
  Graph4->SetLineWidth(4);
  Graph4->SetLineStyle(1);
  Graph4->SetMarkerColor(6);
  Graph4->SetMarkerStyle(21);

  Graph7->SetLineColor(8);
  Graph7->SetLineWidth(4);
  Graph7->SetLineStyle(1);
  Graph7->SetMarkerColor(8);
  Graph7->SetMarkerStyle(21);

  Graph10->SetLineColor(9);
  Graph10->SetLineWidth(4);
  Graph10->SetLineStyle(1);
  Graph10->SetMarkerColor(9);
  Graph10->SetMarkerStyle(21);

  Graph13->SetLineColor(28);
  Graph13->SetLineWidth(4);
  Graph13->SetLineStyle(1);
  Graph13->SetMarkerColor(28);
  Graph13->SetMarkerStyle(21);

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
  Graph4->SetTitle("");
  Graph4->GetXaxis()->SetTitle("#tau#rightarrow3#mu trigger efficiency");
  Graph4->GetYaxis()->SetTitle("MinBias background rejection");
  Graph4->GetXaxis()->SetLabelFont(42);
  Graph4->GetYaxis()->SetLabelFont(42);
  Graph4->GetYaxis()->SetTitleFont(42);
  Graph4->GetXaxis()->SetTitleFont(42);     
  Graph4->GetYaxis()->SetTitleOffset(1.5);
  Graph4->GetYaxis()->SetTitleSize(0.055);
  Graph4->GetYaxis()->SetLabelSize(0.040);
  Graph4->GetXaxis()->SetTitleSize(0.045);
  Graph4->GetXaxis()->SetLabelSize(0.040);
  Graph4->Draw("ALP");
  Graph4->GetXaxis()->SetRangeUser(0,0.4);
  Graph4->GetYaxis()->SetRangeUser(0.99967,1);
  //Graph5->Draw("SAMELP");
  //Graph6->Draw("SAMELP");
  Graph7->Draw("SAMElp");
  //Graph8->Draw("SAME");
  //Graph9->Draw("SAME");
  Graph10->Draw("SAMElP");
  ///Graph11->Draw("SAMELP");
  //Graph12->Draw("SAMELP");
  Graph13->Draw("SAMElp");
  //Graph14->Draw("SAME");
  //Graph15->Draw("SAME");

  TLine *line=new TLine(0, 1-(10/(2760*11.246)), 0.4, 1-(10/(2760*11.246)));
  line->SetLineWidth(7); line->SetLineStyle(1); line->SetLineColor(kOrange);
  line->Draw("SAME");
  c_1->Update();  
  TLegend *leg_1 = new TLegend(0.50,0.70,0.85,0.85,NULL,"brNDC");
  leg_1->SetNColumns(1);leg_1->SetBorderSize(0);
  leg_1->SetTextSize(0.03);leg_1->SetTextFont(42);
  leg_1->SetLineColor(0);leg_1->SetLineWidth(10);
  leg_1->SetFillColor(kWhite);leg_1->SetFillStyle(1001);
  leg_1->AddEntry(Graph4,"L1T_2TkMuStub1Tk", "lp"); 
  //leg_1->AddEntry(Graph5,"L1T_2TkMuStub1Tk_wm0p5_dR", "lp");
  //leg_1->AddEntry(Graph6,"L1T_2TkMuStub1Tk_wm0p7_dR", "lp");
  leg_1->AddEntry(Graph7,"L1T_2TkMu1Tk", "lp"); 
  //leg_1->AddEntry(Graph8,"L1T_2TkMu1Tk_wm0p5_dR", "l");
  //leg_1->AddEntry(Graph9,"L1T_2TkMu1Tk_wm0p7_dR", "l");
  leg_1->AddEntry(Graph10,"L1T_2TkMuStub1TkMu", "lp"); 
  //leg_1->AddEntry(Graph11,"L1T_2TkMuStub1TkMu_wm0p5_dR", "lp");
  //leg_1->AddEntry(Graph12,"L1T_2TkMuStub1TkMu_wm0p7_dR", "lp");
  leg_1->AddEntry(Graph13,"L1T_2TkMu1TkMuStub", "lp"); 
  //leg_1->AddEntry(Graph14,"L1T_2TkMu1TkMuStub_wm0p5_dR", "l");
  //leg_1->AddEntry(Graph15,"L1T_2TkMu1TkMuStub_wm0p7_dR", "l");
  leg_1->AddEntry(line,"Bkg reduction at 10 kHz", "l");
  leg_1->Draw();
  lumi_13TeV = ""; 
  CMS_lumi( p_2, iPeriod, iPos );  
  system ("mkdir myplots"); 
  system ("mkdir myplots/ROCcurves"); 
  char* output;
  if (ext==0)
     output = Form("myplots/ROCcurves/Tau3MuTriggersROC_mix_%i.png",pu) ;
  else
     output = Form("myplots/ROCcurves/Tau3MuTriggersROC_mix_%i_extended.png",pu) ;
  c_1->SaveAs(output);


}

