//  Data unfolding using Singular Value Decomposition 
// 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSVDUnfold example                                                   //
//                                                                      //
// Data unfolding using Singular Value Decomposition (hep-ph/9509307)   //
// Authors: Kerstin Tackmann, Andreas Hoecker, Heiko Lacker             //
// Modified by : Fahmi Maulida                                          //
// Example distribution and smearing model from Tim Adye (RAL)          //
//                                                                      //
// 12-October-2015                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TSVDUnfold.h"
#endif

void UnfoldTSVDptttbar(TString varToUnfold="ptttbar",int kterm=2,int nbinsunfolded=10,int nbinsreconstructed=10, float xmin=0,float xmax=100) 
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  //variables: generated, reconstructed, event weight                                                                 
  Float_t gen_var, rec_var, weight;


  // --- Data/MC ----------------------------------

  // The MC input

  //readout signal expectations and fill response matrix                                                              
  TChain *signal  = new TChain("dataAnalyzer/ue");
  signal->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_TTJets_MSDecays_172v5_9.root");
  signal->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  signal->SetBranchAddress("gen_"+varToUnfold, &gen_var);
  signal->SetBranchAddress("weight", &weight);

  // Histogram for MC                                                                                                 
  TH1D *hTrue = new TH1D("true", "true;"+varToUnfold,     nbinsunfolded, xmin,xmax);
  hTrue->Sumw2();
  TH1D *hMeas = new TH1D("meas", "measured;"+varToUnfold, nbinsunfolded, xmin, xmax);
  hMeas->Sumw2();
  TH2D *DetRes = new TH2D ("DetRes", "covariance matrix", nbinsunfolded, xmin, xmax, nbinsunfolded, xmin, xmax);
  //DetRes->Sum2();


  // Load for Response Matrix                                                                                         
  for (Int_t i= 0; i<signal->GetEntries(); i++) {
    signal->GetEntry(i);
    DetRes->Fill(rec_var, gen_var, weight); //Detector response
    hMeas->Fill(rec_var, weight); //MC Reco
    hTrue->Fill(gen_var,weight); //MC Truth
  }

  // ------Data Input-------------------------------
  TChain *ch_data = new TChain("dataAnalyzer/ue");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012A_0.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012A_1.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012B_0.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012B_1.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012B_2.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012B_3.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012B_4.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_0.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_1.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_2.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_3.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_4.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012C_5.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_0.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_1.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_2.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_3.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_4.root");
  ch_data->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/data/Data8TeV_MuEG2012D_5.root");
  ch_data->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  ch_data->SetBranchAddress("weight", &weight);

  // Histogram for Data                                                                                               
  TH1D *hData  = new TH1D ("Data", "data;"+varToUnfold, nbinsreconstructed, xmin,xmax);
  hData->Sumw2();
  TH1D *hDataTrue  = new TH1D ("Data Truth", "dataTruth;"+varToUnfold, nbinsreconstructed, xmin,xmax);
  hDataTrue->Sumw2();
  TH2D *statcov = new TH2D ("statcov", "covariance matrix", nbinsunfolded, xmin, xmax, nbinsunfolded, xmin, xmax);
   

  for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
    ch_data->GetEntry(i);
    hData->Fill(rec_var); //Data Reco                                                                  
    hDataTrue->Fill(gen_var); // Data truth
  }


  // With BackGround Substract
  //----------------------------------------------Background _WW----------------------------------------------
   
  TChain *bkg_WW = new TChain("dataAnalyzer/ue");
  bkg_WW->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_WW_0.root");
   
  TFile* file=new TFile("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_WW_0.root", "update");
  TVectorD *constVals = (TVectorD *)file->Get("constVals");
  float norigEvents = (*constVals)[0];
  float xsec = (*constVals)[1];
  float lumi = 19701.0;
  float xsecWeight = lumi*xsec/norigEvents;
 
  bkg_WW->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  bkg_WW->SetBranchAddress("weight", &weight);
 
  //histogram for Background
  TH1D * hBkg_WW   = new TH1D ("Bkg_WW",  "bkg_WW;"+varToUnfold,  nbinsreconstructed, xmin,xmax);
  hBkg_WW->Sumw2();
  for (Int_t i= 0; i<bkg_WW->GetEntries(); i++) {
    bkg_WW->GetEntry(i);
    float eventWeight(xsecWeight*weight);
    hBkg_WW->Fill(rec_var, eventWeight);
 
     
  }
  //----------------------------------------------Background _WZ----------------------------------------------                                  
 
  TChain *bkg_WZ = new TChain("dataAnalyzer/ue");
  bkg_WZ->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_WZ_0.root");
  TFile* file=new TFile("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_WZ_0.root", "update");
 
  TVectorD *constVals = (TVectorD *)file->Get("constVals");
  float norigEvents = (*constVals)[0];
  float xsec = (*constVals)[1];
  float lumi = 19701.0;
  float xsecWeight = lumi*xsec/norigEvents;
 
  bkg_WZ->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  bkg_WZ->SetBranchAddress("weight", &weight);
 
  //histogram for Background                                                                                                                    
  TH1D * hBkg_WZ   = new TH1D ("Bkg_WZ",  "bkg_WZ;"+varToUnfold,  nbinsreconstructed, xmin,xmax);
  hBkg_WZ->Sumw2();
  for (Int_t i= 0; i<bkg_WZ->GetEntries(); i++) {
    bkg_WZ->GetEntry(i);
    float eventWeight(xsecWeight*weight);
    hBkg_WZ->Fill(rec_var,eventWeight);
 
  }
 
  //----------------------------------------------Background _ZZ----------------------------------------------                                 
 
  TChain *bkg_ZZ = new TChain("dataAnalyzer/ue");
  bkg_ZZ->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_ZZ_0.root");
  TFile* file=new TFile("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_ZZ_0.root", "update");
  TVectorD *constVals = (TVectorD *)file->Get("constVals");
  float norigEvents = (*constVals)[0];
  float xsec = (*constVals)[1];
  float lumi = 19701.0;
  float xsecWeight = lumi*xsec/norigEvents;
 
 
  bkg_ZZ->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  bkg_ZZ->SetBranchAddress("weight", &weight);
 
  //histogram for Background                                                                                                                    
  TH1D * hBkg_ZZ   = new TH1D ("Bkg_ZZ",  "bkg_ZZ;"+varToUnfold,  nbinsreconstructed, xmin,xmax);
  hBkg_ZZ->Sumw2();
  for (Int_t i= 0; i<bkg_ZZ->GetEntries(); i++) {
    bkg_ZZ->GetEntry(i);
    float eventWeight( xsecWeight*weight );
    hBkg_ZZ->Fill(rec_var,eventWeight);
 
  }
 
  //----------------------------------------------Background _SingleT_tW----------------------------------------------                                                                                                                     
 
  TChain *bkg_singleT_tW = new TChain("dataAnalyzer/ue");
  bkg_singleT_tW->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_SingleT_tW.root");
 
  TFile* file=new TFile("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_SingleT_tW.root", "update");
  TVectorD *constVals = (TVectorD *)file->Get("constVals");
  float norigEvents = (*constVals)[0];
  float xsec = (*constVals)[1];
  float lumi = 19701.0;
  float xsecWeight = lumi*xsec/norigEvents;
 
  bkg_singleT_tW->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  bkg_singleT_tW->SetBranchAddress("weight", &weight);
 
  //histogram for Background                                                                                    
  TH1D * hBkg_singleT_tW   = new TH1D ("Bkg_singleT_TW",  "bkg_singleT_tW;"+varToUnfold,  nbinsreconstructed, xmin,xmax);
  hBkg_singleT_tW->Sumw2();
  for (Int_t i= 0; i<bkg_singleT_tW->GetEntries(); i++) {
    bkg_singleT_tW->GetEntry(i);
 
    float eventWeight( xsecWeight*weight );
    hBkg_singleT_tW->Fill(rec_var,eventWeight);
 
 
  }
 
  //----------------------------------------------Background _SingleTbar_tW----------------------------------------------                                                                                    
  TChain *bkg_singleTbar_tW = new TChain("dataAnalyzer/ue");
  bkg_singleTbar_tW->Add("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_SingleTbar_tW.root");
 
  TFile* file=new TFile("/afs/cern.ch/user/f/fmaulida/CMSSW_7_2_0_pre3/src/Unfolding/RooUnfold-1.1.1/MC/MC8TeV_SingleTbar_tW.root", "update");
  TVectorD *constVals = (TVectorD *)file->Get("constVals");
  float norigEvents = (*constVals)[0];
  float xsec = (*constVals)[1];
  float lumi = 19701.0;
  float xsecWeight = lumi*xsec/norigEvents;
 
  bkg_singleTbar_tW->SetBranchAddress("rec_"+varToUnfold, &rec_var);
  bkg_singleTbar_tW->SetBranchAddress("weight", &weight);
 
  //histogram for Background         
  TH1D * hBkg_singleTbar_tW   = new TH1D ("Bkg_singleTbar_tW",  "bkg_singleTbar_tW;"+varToUnfold,  nbinsreconstructed, xmin,xmax);
  hBkg_singleTbar_tW->Sumw2();
  for (Int_t i= 0; i<bkg_singleTbar_tW->GetEntries(); i++) {
    bkg_singleTbar_tW->GetEntry(i);
 
    float eventWeight( xsecWeight*weight );
    hBkg_singleTbar_tW->Fill(rec_var,eventWeight);
  }
 
  //---------------------------------------Data with background substract------------------------------------------
  //subtract the background contribution from the data
  hData->Add(hBkg_WW,-1);
  hData->Add(hBkg_WZ,-1);
  hData->Add(hBkg_ZZ,-1);
  hData->Add(hBkg_singleT_tW,-1);
  hData->Add(hBkg_singleTbar_tW,-1);

  // Fill the data covariance matrix
  for (int i=1; i<=hData->GetNbinsX(); i++) {
    statcov->SetBinContent(i,i,hData->GetBinError(i)*hData->GetBinError(i)); 
  }

  // --- Here starts the actual unfolding -------------------------

  // Create TSVDUnfold object and initialise
  TSVDUnfold *tsvdunf = new TSVDUnfold(hData, statcov, hMeas, hTrue, DetRes );

  // It is possible to normalise unfolded spectrum to unit area
  //  tsvdunf->SetNormalize(kFALSE); // no normalisation here

  // Perform the unfolding with regularisation parameter kreg = 13
  // - the larger kreg, the finer grained the unfolding, but the more fluctuations occur
  // - the smaller kreg, the stronger is the regularisation and the bias
  TH1D* unfres = tsvdunf->Unfold( kterm );

  // Get the distribution of the d to cross check the regularization
  // - choose kreg to be the point where |d_i| stop being statistically significantly >>1
  TH1D* ddist = tsvdunf->GetD();

  // Get the distribution of the singular values
  TH1D* svdist = tsvdunf->GetSV();

  // Compute the error matrix for the unfolded spectrum using toy MC
  // using the measured covariance matrix as input to generate the toys
  // 100 toys should usually be enough
  // The same method can be used for different covariance matrices separately.
  TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix( statcov, 100 );   

  // Now compute the error matrix on the unfolded distribution originating
  // from the finite detector matrix statistics
  TH2D* uadetcov = tsvdunf->GetAdetCovMatrix( 100 );   

  // Sum up the two (they are uncorrelated)
  ustatcov->Add( uadetcov );

  //Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics. 
  TH2D* utaucov = tsvdunf->GetXtau();
  utaucov->Add( uadetcov );

  //Get the computed inverse of the covariance matrix
  TH2D* uinvcov = tsvdunf->GetXinv();

   
  // --- Only plotting stuff below ------------------------------
   
  
  for (int i=1; i<=unfres->GetNbinsX(); i++) {
    unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i,i)));
  }

  // Renormalize just to be able to plot on the same scale
  hTrue->Scale(0.7*hDataTrue->Integral()/hTrue->Integral());

  TLegend *leg = new TLegend(0.58,0.68,0.99,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(unfres,"Unfolded Data","p");
  //  leg->AddEntry(hDataTrue,"True Data","l");
  leg->AddEntry(hData,"Reconstructed Data","l");
  leg->AddEntry(hTrue,"True MC","l");
  leg->AddEntry(hBkg_WW,"WW","l");
  leg->AddEntry(hBkg_WZ,"WZ","l");
  leg->AddEntry(hBkg_ZZ,"ZZ","l");
  leg->AddEntry(hBkg_singleT_tW,"Single Top","l");
  leg->AddEntry(hBkg_singleTbar_tW,"Single Tbar","l");

  TCanvas *c1 = new TCanvas( "c1","ptttbar Unfolding with TSVDUnfold ", 900, 800 );

  // --- Style settings -----------------------------------------
  Int_t c_Canvas    = TColor::GetColor( "#f0f0f0" );
  Int_t c_FrameFill = TColor::GetColor( "#fffffd" );
  Int_t c_TitleBox  = TColor::GetColor( "#6D7B8D" );
  Int_t c_TitleText = TColor::GetColor( "#FFFFFF" );

  c1->SetFrameFillColor( c_FrameFill );
  c1->SetFillColor     ( c_Canvas    );
  c1->Divide(1,2);
  TVirtualPad * c11 = c1->cd(1);
  c11->SetFrameFillColor( c_FrameFill );
  c11->SetFillColor     ( c_Canvas    );

  gStyle->SetTitleFillColor( c_TitleBox  );
  gStyle->SetTitleTextColor( c_TitleText );
  gStyle->SetTitleBorderSize( 1 );
  gStyle->SetTitleH( 0.052 );
  gStyle->SetTitleX( c1->GetLeftMargin() );
  gStyle->SetTitleY( 1 - c1->GetTopMargin() + gStyle->GetTitleH() );
  gStyle->SetTitleW( 1 - c1->GetLeftMargin() - c1->GetRightMargin() );

  TH1D* frame = new TH1D( *unfres );
  frame->SetTitle("ptttbar Unfolding with with TSVDUnfold" );
  frame->GetXaxis()->SetTitle( "x variable" );
  frame->GetYaxis()->SetTitle( "Events" );
  frame->GetXaxis()->SetTitleOffset( 1.25 );
  frame->GetYaxis()->SetTitleOffset( 1.29 );
  frame->Draw();


  //--------------------------------------                                                                       
  hBkg_WW->SetLineStyle(2);
  hBkg_WW->SetLineColor(2);
  hBkg_WW->SetLineWidth(2);
  //-------------------------------------                                                                        
  hBkg_WZ->SetLineStyle(2);
  hBkg_WZ->SetLineColor(9);
  hBkg_WZ->SetLineWidth(2);
  //------------------------------------                                                                         
  hBkg_ZZ->SetLineStyle(2);
  hBkg_ZZ->SetLineColor(5);
  hBkg_ZZ->SetLineWidth(2);
  //------------------------------------------                                                                   
  hBkg_singleT_tW->SetLineStyle(2);
  hBkg_singleT_tW->SetLineColor(6);
  hBkg_singleT_tW->SetLineWidth(2);
  //------------------------------------------                                                                   
  hBkg_singleTbar_tW->SetLineStyle(2);
  hBkg_singleTbar_tW->SetLineColor(3);
  hBkg_singleTbar_tW->SetLineWidth(2);


  hData->SetLineStyle(2);
  hData->SetLineColor(4);
  hData->SetLineWidth(2);
  unfres->SetMarkerStyle(20);
  // hDataTrue->SetLineColor(2);
  //hDataTrue->SetLineWidth(2);
  hTrue->SetLineStyle(2);
  hTrue->SetLineColor(8);
  hTrue->SetLineWidth(2);
  // ------------------------------------------------------------

  // add histograms
  unfres->DrawNormalized();
  //hDataTrue->Draw("same");
  hData->DrawNormalized("same");
  hTrue->DrawNormalized("same");
  hBkg_WW->DrawNormalized("same");
  hBkg_WZ->DrawNormalized("same");
  hBkg_ZZ->DrawNormalized("same");
  hBkg_singleT_tW->DrawNormalized("same");
  hBkg_singleTbar_tW->DrawNormalized("same");
  leg->Draw();

  // covariance matrix
  gStyle->SetPalette(1,0);
  TVirtualPad * c12 = c1->cd(2);
  c12->Divide(2,1);
  TVirtualPad * c2 = c12->cd(1);
  c2->SetFrameFillColor( c_FrameFill );
  c2->SetFillColor     ( c_Canvas    );
  c2->SetRightMargin   ( 0.15         );

  TH2D* covframe = new TH2D( *ustatcov );
  covframe->SetTitle( "TSVDUnfold covariance matrix" );
  covframe->GetXaxis()->SetTitle( "x variable" );
  covframe->GetYaxis()->SetTitle( "x variable" );
  covframe->GetXaxis()->SetTitleOffset( 1.25 );
  covframe->GetYaxis()->SetTitleOffset( 1.29 );
  covframe->Draw();

  ustatcov->SetLineWidth( 2 );
  ustatcov->Draw( "colzsame" );

  // distribution of the d quantity
  TVirtualPad * c3 = c12->cd(2);
  c3->SetFrameFillColor( c_FrameFill );
  c3->SetFillColor     ( c_Canvas    );
  c3->SetLogy();

  TLine *line = new TLine( 0.,1.,40.,1. );
  line->SetLineStyle(2);

  TH1D* dframe = new TH1D( *ddist );
  dframe->SetTitle( "TSVDUnfold |d_{i}|" );
  dframe->GetXaxis()->SetTitle( "i" );
  dframe->GetYaxis()->SetTitle( "|d_{i}|" );
  dframe->GetXaxis()->SetTitleOffset( 1.25 );
  dframe->GetYaxis()->SetTitleOffset( 1.29 );
  dframe->SetMinimum( 0.001 );
  dframe->Draw();

  ddist->SetLineWidth( 2 );
  ddist->Draw( "same" );
  line->Draw();
}
