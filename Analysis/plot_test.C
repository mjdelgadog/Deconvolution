#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

int plot_test(const char* file0, const char* file1, int iev = 1, int ich = 456) {
  TFile* file[2] = {0}; 
  file[0] = new TFile(file0); 
  file[1] = new TFile(file1); 

  TH1D* hWaveform[2] = {0}; 
  TH1D* hDeconv  [2] = {0};
  TH1D* hSDFilter[2] = {0}; 

  Color_t col[2] = {kBlue+1, kRed+1}; 

  for(int i=0; i<2; i++) {
    printf("Getting waveform...\n"); 
    hWaveform[i] = (TH1D*)file[i]->Get(
        Form("opdecoanav/event_%i_opchannel_%i_decowaveform_1", 
          iev, ich)); 
    printf("Getting deconvolution...\n"); 
    hDeconv[i] = (TH1D*)file[i]->Get(
        Form("opdecoana/event_%i_opchannel_%i_decowaveform_1", 
          iev, ich)); 
    printf("Getting filter spectral density...\n"); 
    hSDFilter[i] = (TH1D*)file[i]->Get(
        Form("opdecoanaG/event_%i_opchannel_%i_decowaveform_1", 
          iev, ich)); 

    hWaveform[i]->SetLineColor(col[i]); 
    hWaveform[i]->SetLineWidth(2); 

    hDeconv[i]->SetLineColor(col[i]); 
    hDeconv[i]->SetLineWidth(2); 

    hSDFilter[i]->SetLineColor(col[i]); 
    hSDFilter[i]->SetLineWidth(2); 
  }

  TCanvas* cPlot = new TCanvas("cPlot", "ev plot", 0, 0, 900, 1000); 
  TPad* pWaveform = new TPad("pWaveform", "Waveform", 0, 0.6, 0.5, 1); 
  TPad* pDeconv   = new TPad("pDeconv", "Deconvolution", 0, 0.2, 0.5, 0.6); 
  TPad* pResidual = new TPad("pResidual", "Residuals", 0, 0, 0.5, 0.2); 
  TPad* pSD       = new TPad("pSD", "Spectral Density", 0.5, 0., 1.0, 1.0); 


  pWaveform->Draw(); pWaveform->cd(); 
  hWaveform[0]->Draw("hist"); 
  hWaveform[1]->Draw("hist same"); 

  cPlot->cd(); pDeconv->Draw(); pDeconv->cd();  
  hDeconv[0]->Draw("hist"); 
  hDeconv[1]->Draw("hist same"); 

  cPlot->cd(); pResidual->Draw(); pResidual->cd();  
  TH1D* hRes = (TH1D*) hDeconv[0]->Clone("Residuals"); 
  hRes->Add(hDeconv[1], -1); 
  hRes->Draw(); 

  cPlot->cd(); pSD->Draw(); pSD->cd();  
  hSDFilter[0]->Draw("hist"); 
  hSDFilter[1]->Draw("hist same"); 


  return 1; 
}
