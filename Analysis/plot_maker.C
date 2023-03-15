/**
 * @author      : dguffant (dguffant@dunegpvm11.fnal.gov)
 * @file        : plot_maker
 * @created     : lunedì mar 13, 2023 05:29:55 CDT
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TPRegexp.h"
#include "TTimer.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"

std::vector<float> compute_derivative(const std::vector<float>& v, const float dx = 1.0) {
  std::vector<float> vout(v.size(), 0.);
  vout.at(0) = (v.at(1) - v.at(0)) / dx; 
  vout.at(1) = (v.at(2) - v.at(0)) / (2*dx); 
  vout.at(v.size()-2) = (v.at(v.size()-1) -  v.at(v.size()-3)) / (2*dx); 
  vout.back () = (v.back() - v.at(v.size() - 2))/dx; 

  float den = 8*dx; 
  for (int i=2; i<v.size() - 2; i++) {
    float num = -v.at(i+2) +8*v.at(i+1) -8*v.at(i-1) + v.at(i-2); 
    vout.at(i) = num / den; 
  }
  return vout; 
}

std::vector<float> compute_derivative(const TH1D* h, const float dx = 1) {
  std::vector<float> v(h->GetNbinsX(), 0.); 
  for (int i=1; i<=h->GetNbinsX(); i++) v.at(i-1) = h->GetBinContent(i); 
  auto vout = compute_derivative(v, dx); 
  return vout; 
}

TPolyMarker* find_peaks(const TH1D* hin, const double thrs) {
  auto der = compute_derivative(hin, 1.0);
  auto dder= compute_derivative(der, 1.0);

  std::vector<double> xpos;
  std::vector<double> xval; 
  TPolyMarker* pm = nullptr; 

  for (int i=1; i<der.size(); i++) {
    if (der.at(i) * der.at(i-1) <= 0) {
      if (dder.at(i) > 0) continue;

      if (hin->GetBinContent(i+1) >= thrs) {
        xpos.push_back(hin->GetBinCenter(i+1)); 
        xval.push_back(hin->GetBinContent(i+1)); 
      }
    }
  }

  if (!xpos.empty()) {
    pm = new TPolyMarker(xpos.size(), &xpos[0], &xval[0]); 
    pm->SetMarkerStyle(23); 
    pm->SetMarkerColor(kOrange+7); 
  }

  return pm; 
}

float find_pe_by_charge(const TH1* hin, const TPolyMarker* peaks, std::vector<TH1*>* hintegral = 0) {
  int   i_lower_limit = 0;
  float x_lower_limit = 0; 
  int   i_upper_limit = 0; 
  float x_upper_limit = 0; 

  float charge = 0.;
  for (int j=0; j<peaks->GetN(); j++) {
    float x_peak = peaks->GetX()[j]; 
    // ignore if peak position is within the already explored waveform interval
    if (x_peak < x_upper_limit) continue;

    int i_peak = hin->GetXaxis()->FindBin(x_peak); 
    // compute integral moving left
    int i = i_peak;
    while (hin->GetBinContent(i) >= 0 && i>1) {
      charge += hin->GetBinContent(i); i--;
    }
    i_lower_limit = i; 
    x_lower_limit = hin->GetXaxis()->GetBinCenter(i); 

    i = i_peak; 
    while (hin->GetBinContent(i) >= 0 && i<hin->GetNbinsX()) {
      charge += hin->GetBinContent(i); i++;
    }
    i_upper_limit = i;
    x_upper_limit = hin->GetXaxis()->GetBinCenter(i); 

    if (hintegral) {
      TH1* h = (TH1*)hin->Clone(); 
      h->SetNameTitle(Form("%s_pk_%i", hin->GetName(), j), Form("%s - pulse %i", hin->GetTitle(), j)); 
      h->Reset(); 
      for (int k=1; k<h->GetNbinsX(); k++) {
        if (k>=i_lower_limit && k<=i_upper_limit) {
          h->SetBinContent(k, hin->GetBinContent(k)); 
        }
      }

      hintegral->push_back(h); 
    }
  }

  return charge; 
}

void plot_maker(const char* raw_file, const char* wdec_file_pf, const char* wdec_file_npf = "") {
  int iBLPretrigger = 30; 
  TFile* raw = new TFile(raw_file); 
  TFile* wdec_pf = new TFile(wdec_file_pf); 
  TFile* wdec_npf = new TFile(wdec_file_npf); 

  TTree* traw = (TTree*)raw->Get("opdigi/PhotonData"); 
  std::vector<int> *raw_pe = 0; 
  std::vector<double> *raw_time = 0; 
  traw->SetBranchAddress("photon_opCh", &raw_pe); 
  traw->SetBranchAddress("photon_pulse", &raw_time); 


  TH2D* h2ChTime = new TH2D("h2ChTime", "raw pe;Ch number;Time [#mu s]", 
      500, 0.5, 500.5, 1000, 0, 16); 

  traw->GetEntry(0); 
  assert(raw_pe->size() == raw_time->size()); 

  for (size_t i=0; i<raw_pe->size(); i++) {
    h2ChTime->Fill(raw_pe->at(i), raw_time->at(i)); 
  }

  TCanvas* raw_c = new TCanvas("raw_c", "raw_c", 0, 0, 800, 800);
  raw_c->cd(); 
  h2ChTime->Draw("colz"); 

  TDirectory* opdecoanapf = (TDirectory*)wdec_pf->Get("opdecoana"); 
  //TDirectory* opdecoananpf = (TDirectory*)wdec_npf->Get("opdecoana"); 
  TPRegexp rxh("event_\\d+_opchannel_(\\d+)_decowaveform_(\\d+)"); 

  TFile* analysis_output = new TFile("analysis_output.root", "recreate"); 
  TTree* t = new TTree("t", "deco analysis tree"); 
  double true_pe = 0; 
  double wdec_pe = 0; 
  double wpulse_pe = 0; 

  t->Branch("true_pe", &true_pe);
  t->Branch("wdec_pe", &wdec_pe);
  t->Branch("wpulse_pe", &wpulse_pe);

  TCanvas* ctest = new TCanvas("ctest", "test", 0, 0, 800, 700); 
  bool vis = true; 

  auto shift_hist = [](const TH1* h, const int ishift) {
    TH1D* hh = (TH1D*)h->Clone(); 
    hh->Reset(); 
    for (int i=ishift+1; i<=h->GetNbinsX()-ishift; i++) {
      //printf("origin: %i -> target = %i\n", i-ishift, i); 
      hh->SetBinContent(i, h->GetBinContent(i-ishift)); 
    }
    return hh; 
  };

  TTimer* timer = new TTimer("gSystem->ProcessEvents();", 500, false); 
  TH1D* hBaselineDec = new TH1D("hBaselineDec", Form("Baseline dec: [0-%i] samples;Waveform values [a.u];Entries", iBLPretrigger), 
      200, -10, +10); 
  TH1D* hBaselineRaw = new TH1D("hBaselineRaw", "Baseline raw: full sample;Waveform values [a.u];Entries", 
      200, -10, +10); 

  for (const auto &dec_ : *opdecoanapf->GetListOfKeys()) {
    TH1D* hwdec_pf = (TH1D*)wdec_pf ->Get(Form("opdecoana/%s", dec_->GetName()));
    TH1D* hwdec_npf= (TH1D*)wdec_npf->Get(Form("opdecoana/%s", dec_->GetName()));

    auto match = rxh.MatchS(hwdec_pf->GetName()); 
    
    TObjString* osCh = (TObjString*)match->At(1); 
    TObjString* osWv = (TObjString*)match->At(2); 
    int iCh = osCh->GetString().Atoi(); 
    int iWv = osWv->GetString().Atoi(); 
    int iWv_raw = iWv - 1; if (iWv == 3) iWv_raw = 1;  
    TString hraw_name = Form("opdigiana/event_1_opchannel_%i_waveform_%i", iCh, iWv_raw); 
    printf("event 1: opchannel %i - wv %i: %s\n", iCh, iWv, hraw_name.Data()); 

    TH1D* hraw = (TH1D*)raw->Get(hraw_name); 

    double tmin = hraw->GetXaxis()->GetXmin(); 
    double tmax = hraw->GetXaxis()->GetXmax(); 
    TH1D* htrue = h2ChTime->ProjectionY("hraw", iCh, iCh, ""); 
    int itimemin = htrue->GetXaxis()->FindBin(tmin);  
    int itimemax = htrue->GetXaxis()->FindBin(tmax);  
    true_pe = htrue->Integral(itimemin, itimemax); 
    wdec_pe = hwdec_pf->Integral(); 
    wpulse_pe = 0.; 

    double spt_thrs = TMath::Min(0.9, 0.15 / hwdec_pf->GetMaximum()); 
    auto peaks = find_peaks(hwdec_pf, 0.15); 

    if (true_pe>500) {
      for (int i=1; i<=iBLPretrigger; i++) hBaselineDec->Fill( hwdec_pf->GetBinContent(i) ); 
    } else if (true_pe == 0) {
      for (int i=1; i<=hraw->GetNbinsX(); i++) hBaselineRaw->Fill( hraw->GetBinContent(i) - 1500 ); 
    }
    auto hintegral = new std::vector<TH1*>;
    if (peaks) {
      printf("%i peaks found - thrs: %g\n", peaks->GetN(), spt_thrs); 
      wpulse_pe = find_pe_by_charge(hwdec_pf, peaks/*, hintegral*/);
    }
    else printf("No peaks found\n"); 

    t->Fill(); 
    bool do_break = false; 

    if (vis) {
      timer->TurnOn(); 
      timer->Reset(); 

      TH1D* htrue_shift = shift_hist(htrue, 96); 
      htrue_shift->SetFillColor(kGray+1); htrue_shift->SetLineWidth(0); 
      ctest->Clear(); 
      ctest->Divide(1, 2); 
      ctest->cd(1); hraw->Draw(); gPad->Modified(); gPad->Update(); 
      ctest->cd(2); 
      hwdec_npf->SetLineColorAlpha(kAzure, 0.2); 
      hwdec_npf->Draw("hist"); 
      //hwdec_pf->Draw("axis"); 
      htrue_shift->Draw("hist same");  
      for (auto &h : *hintegral) {
        h->SetFillColor(kAzure-4); 
        h->SetLineWidth(0); 
        h->Draw("hist same"); 
      }
      hwdec_pf->SetLineColor(kBlack); hwdec_pf->SetLineWidth(2); 
      hwdec_pf->Draw("hist same");
      if (peaks) peaks->Draw("same"); 

      gPad->Modified(); gPad->Update(); 
      printf("true pe: %g, wdec pe: %g, pulse pe: %g\n", true_pe, wdec_pe, wpulse_pe); 
      printf("press 'q' to close loop, any other key to continue\n"); 
      char input_c; 
      std::cin >> input_c; 
      if (input_c == 'q') do_break = true; 
      timer->TurnOff(); 
    }
    
    if (hintegral) delete hintegral; 
    delete htrue; 

    if (do_break) break; 
  }
  TCanvas* cBaseline = new TCanvas("cBaseline"); 
  cBaseline->cd(); 
  hBaselineDec->Draw(); 
  hBaselineRaw->Draw("hist same"); 

  //TH1D* hdiff = new TH1D("hdiff", "rec-true;(pe reco - pe truth)/pe truth;Entries",
      //100, -100, +100); 
  //t->Draw("(wdec_pe-true_pe)/true_pe*100", "true_pe>0", ""); 
  //ctest->Clear(); 
  //hdiff->DrawClone("hist"); 
  //ctest->Update(); 

  t->Write(); 
  //analysis_output->Close(); 

  return; 
}



