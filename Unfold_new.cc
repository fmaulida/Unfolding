#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
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
//#include "TUnfoldSys.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TSVDUnfold.h"
#endif
#include "TUnfoldDensity.h"
#include "TChain.h"

void Unfold_new()
{
  TString varOutput, varToUnfold, varReg, title, Reg;
  Int_t reg = 4;
  TString input = "/afs/cern.ch/work/f/fmaulida/13TeV/CMSSW_8_0_5_patch1/src/Unfolding/newlooks/binning/";
  TString output = "/afs/cern.ch/work/f/fmaulida/13TeV/CMSSW_8_0_5_patch1/src/Unfolding/newlooks/UE_16_4/";

  cout<<"varToUnfold:"<<endl;
  cin>>varToUnfold;

  cout<<"varReg:"<<endl;
  cin>>varReg;

  cout<<"Output:"<<endl;
  cin>>varOutput;

  if(varToUnfold=="avgptflux") title = "<#sum p_{T}> ";
  if(varToUnfold=="ptflux")    title = "#sum p_{T}";
  if(varToUnfold=="nch")       title = "N_{ch}";

  if(varReg=="inc")            Reg = "(Inclusive)";
  if(varReg=="tran")           Reg = "(Transverse)";
  if(varReg=="away")           Reg = "(Away)";
  if(varReg=="tow")            Reg = "(Toward)";

  //========Input_MC_binning===========
  TH1D *hMC_rec;
  TH1D *hnewMC_rec;
  TH1D *hMC_gen;
  TH1D *hnewMC_gen;
  TH2D *hMC_genrec;
  TH2D *hnewMC_genrec;
 
  TFile *MC_file = new TFile(input+"out_MC_binning.root");
  
  {
    hMC_gen=(TH1D*)MC_file->Get("Gen_"+varToUnfold+"_"+varReg);
    hnewMC_gen=(TH1D*)MC_file->Get("newGen_"+varToUnfold+"_"+varReg);
    hMC_rec=(TH1D*)MC_file->Get("Rec_"+varToUnfold+"_"+varReg);
    hnewMC_rec=(TH1D*)MC_file->Get("newRec_"+varToUnfold+"_"+varReg);
    hMC_genrec=(TH2D*)MC_file->Get("GenRec_"+varToUnfold+"_"+varReg+"_"+varReg);
    hnewMC_genrec=(TH2D*)MC_file->Get("newGenRec_"+varToUnfold+"_"+varReg+"_"+varReg);
  }


  //=========Data=============================
  TH1D *hData_rec;
  TH1D *hnewData_rec;
  
  TFile *data_file = new TFile(input+"out_data_binning.root");
  
  {
    hData_rec=(TH1D*)data_file->Get("Rec_data_"+varToUnfold+"_"+varReg);
    hnewData_rec=(TH1D*)data_file->Get("new_Rec_data"+varToUnfold+"_"+varReg);      
  }
  
  //=======bkg=================
  TH1D *hBkg_WW;
  TH1D *hBkg_WZ;
  TH1D *hBkg_ZZ;
  TH1D *hBkg_singleT_tW;
  TH1D *hBkg_singleTbar_tW;

  TFile *bkg_file = new TFile(input+"out_bkg_binning.root");

  {
    hBkg_WW=(TH1D*)bkg_file->Get("Bkg_WW_"+varToUnfold+"_"+varReg);
    hBkg_WZ=(TH1D*)bkg_file->Get("Bkg_WZ_"+varToUnfold+"_"+varReg);
    hBkg_ZZ=(TH1D*)bkg_file->Get("Bkg_ZZ_"+varToUnfold+"_"+varReg);
    hBkg_singleT_tW=(TH1D*)bkg_file->Get("Bkg_singleT_tW_"+varToUnfold+"_"+varReg);
    hBkg_singleTbar_tW=(TH1D*)bkg_file->Get("Bkg_singleTbar_tW_"+varToUnfold+"_"+varReg);
  }

  //========Subtract the data===============
  TH1D *hData_subtracted;
  TH1D *hnewData_subtracted;
  
  {  
    hData_subtracted=(TH1D*)data_file->Get("Data_sub_"+varToUnfold+"_"+varReg);
    hnewData_subtracted=(TH1D*)data_file->Get("newData_sub_"+varToUnfold+"_"+varReg);
  }
  
  //=======Unfold Data subtracted=======================
  TH1 *hData_unfolded;
  TH1 *hData_sub_unfolded;
  TH1 *hRec_unfolded;
  TH1 *hData_sub_folded;

  TGraph *rhoAvgGr;
  TGraph *rhoMaxGr;
  Float_t rhoAvg;
  Float_t rhoMax;
  Double_t sub_tauUsed, rec_tauUsed, wi_tauUsed;
  Double_t scaleBias = 1;

  {
    hData_unfolded= (TH1D*)hnewMC_gen->Clone("Rec_Data_unfolded_"+varToUnfold+"_"+varReg);
    hData_sub_unfolded= (TH1D*)hnewMC_gen->Clone("Rec_Data_sub_unfolded_"+varToUnfold+"_"+varReg);
    hRec_unfolded = (TH1D*)hnewMC_gen->Clone("Rec_MC_unfolded_"+varToUnfold+"_"+varReg);  
    hData_sub_folded= (TH1D*)hnewData_subtracted->Clone("Rec_Data_folded_"+varToUnfold+"_"+varReg);
  
    
    TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
    TUnfold::EConstraint constraint = TUnfold::kEConstraintNone;
    TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
    
    //====Data rec without subtracted===============
    TUnfoldDensity wi_unfold(hnewMC_genrec, TUnfold::kHistMapOutputHoriz, regMode, constraint, densityFlags);
    wi_unfold.SetInput(hnewData_rec,scaleBias);
    
    //=== Data rec subtracted ===============
    TUnfoldDensity sub_unfold(hnewMC_genrec, TUnfold::kHistMapOutputHoriz, regMode, constraint, densityFlags);
    sub_unfold.SetInput(hnewData_subtracted,scaleBias);
   
    //===== rec MC ========
    TUnfoldDensity rec_unfold(hnewMC_genrec, TUnfold::kHistMapOutputHoriz, regMode, constraint, densityFlags);
    rec_unfold.SetInput(hnewMC_rec,scaleBias);

        //============ scan tau regulation ================================================================================                                                                       
    Double_t tauMin=3.0e-06;
    Double_t tauMax=2.0e-01;

    //Double_t tauMin=1.e-4;    
    //Double_t tauMax=1.e-1;                                                                
    
    TSpline *logTauX,*logTauY;
    TSpline *sub_rhoScan = 0;
    TSpline *rec_rhoScan = 0;
    TSpline *wi_rhoScan = 0;
    rhoAvgGr=new TGraph();
    rhoMaxGr=new TGraph();
    //Int_t nScan=100;
    Int_t nScan=500;
 
    TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoAvg;
    
    wi_unfold.ScanTau(nScan,tauMin,tauMax,&wi_rhoScan,tauflag);
    int wi_iBest = wi_unfold.ScanTau(nScan,tauMin,tauMax,&wi_rhoScan,tauflag);

    sub_unfold.ScanTau(nScan,tauMin,tauMax,&sub_rhoScan,tauflag);
    int sub_iBest = sub_unfold.ScanTau(nScan,tauMin,tauMax,&sub_rhoScan,tauflag);

    rec_unfold.ScanTau(nScan,tauMin,tauMax,&sub_rhoScan,tauflag);
    int rec_iBest = rec_unfold.ScanTau(nScan,tauMin,tauMax,&sub_rhoScan,tauflag);
  
    wi_tauUsed = wi_unfold.GetTau();
    sub_tauUsed = sub_unfold.GetTau();
    rec_tauUsed = rec_unfold.GetTau();
    
    wi_unfold.DoUnfold(wi_tauUsed,hnewData_rec);
    sub_unfold.DoUnfold(sub_tauUsed,hnewData_subtracted);
    rec_unfold.DoUnfold(rec_tauUsed,hnewMC_rec);

    hData_unfolded=wi_unfold.GetOutput("Rec_Data_unfolded_"+varToUnfold+"_"+varReg);    
    hData_sub_unfolded=sub_unfold.GetOutput("Rec_Data_sub_unfolded_"+varToUnfold+"_"+varReg);    
    hRec_unfolded=rec_unfold.GetOutput("Rec_MC_unfolded_"+varToUnfold+"_"+varReg);    
    
    //====Folded subtracted unfolded data===
    hData_sub_folded=sub_unfold.GetFoldedOutput("Rec_Data_sub_folded_"+varToUnfold+"_"+varReg);    
  
   
    //====Create graphs with one point to visualise best choice of tau
    Double_t wi_t[1],wi_rho[1];
    Double_t sub_t[1],sub_rho[1];//x[1],y[1];
    Double_t rec_t[1],rec_rho[1];//x[1],y[1];
     
    wi_rhoScan->GetKnot(wi_iBest,sub_t[0],wi_rho[0]);
    TGraph *wi_bestRho=new TGraph(1,wi_t,wi_rho);
    Double_t *wi_tAll=new Double_t[nScan],*wi_rhoAll=new Double_t[nScan];
    
    for(Int_t i=0;i<nScan;i++) {
      wi_rhoScan->GetKnot(i,wi_tAll[i],wi_rhoAll[i]);
    }
    TGraph *wi_knots=new TGraph(nScan,wi_tAll,wi_rhoAll);    
    TCanvas *wi_tcanv = new TCanvas();
    wi_tcanv->cd();
    wi_rhoScan->Draw();
    wi_rhoScan->SetTitle(";log_{10}(#tau);average(#rho_{i})");
    wi_rhoScan->SetLineColor(kRed);
    // wi_bestRho->Draw("*");
    wi_tcanv->SaveAs(output+"wi_knots_"+varToUnfold+"_"+varReg+".png");
   
    
    sub_rhoScan->GetKnot(sub_iBest,sub_t[0],sub_rho[0]);
    TGraph *sub_bestRho=new TGraph(1,sub_t,sub_rho);
    Double_t *sub_tAll=new Double_t[nScan],*sub_rhoAll=new Double_t[nScan];

    for(Int_t i=0;i<nScan;i++) {
      sub_rhoScan->GetKnot(i,sub_tAll[i],sub_rhoAll[i]);
    }
    TGraph *sub_knots=new TGraph(nScan,sub_tAll,sub_rhoAll);

    TCanvas *sub_tcanv = new TCanvas();
    sub_tcanv->cd();
    sub_rhoScan->Draw();
    sub_rhoScan->SetTitle(";log_{10}(#tau);average(#rho_{i})");
    sub_rhoScan->SetLineColor(kRed);
    //    sub_bestRho->Draw("*");
    sub_tcanv->SaveAs(output+"sub_knots_"+varToUnfold+"_"+varReg+".png");
    /*

    
    rec_rhoScan->GetKnot(rec_iBest,rec_t[0],rec_rho[0]);
    TGraph *rec_bestRho=new TGraph(1,rec_t,rec_rho);
    Double_t *rec_tAll=new Double_t[nScan],*rec_rhoAll=new Double_t[nScan];
   
    for(Int_t i=0;i<nScan;i++) {
      rec_rhoScan->GetKnot(i,rec_tAll[i],rec_rhoAll[i]);
      //       cout<<i<<"  "<<tAll[i]<<"  "<<rhoAll[i]<<endl;                                                                                                                                 
    }
    
    // cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  TAU,iBest --> "<<unfold.GetTau()<<"  "<<iBest<<endl;                                                                
    TGraph *rec_knots=new TGraph(nScan,rec_tAll,rec_rhoAll);
    
    TCanvas *rec_tcanv = new TCanvas();
    rec_tcanv->cd();
    rec_rhoScan->Draw();
    rec_rhoScan->SetTitle(";log_{10}(#tau);average(#rho_{i})");
    rec_rhoScan->SetLineColor(kRed);
    //rec_bestRho->Draw("*");
    rec_tcanv->SaveAs(output+"rec_knots_"+varToUnfold+"_"+varReg+".png");
    */
  }
      
 
  //====Draw=========================
  //============Rec Data vs data subtracted bkg==========
  gStyle->SetOptStat(0);
  TCanvas *recVsSub_bef;
  TCanvas *recVsSub_aft;
  TCanvas *unfoldedVsGen;
  TCanvas *unfoldedSubVsGen;
  TCanvas *unfoldedVsRecVsGen;
  TCanvas *foldedVsSubData;
  Double_t x1 = 0.2, y1 = 0.75, x2 = 0.4, y2 = 0.8;

  {   //========rec data VS Subtracted data Before Binning opt
    recVsSub_bef= new TCanvas("recVSbkg_bef_"+varToUnfold+"_"+varReg,"rec_Vs_Subtracted_data_"+varToUnfold+"_"+varReg,1000,1000);
    hData_rec->SetLineWidth(2);
    hData_rec->Draw();
    hData_rec->SetTitle("Rec Data Vs Subtracted Data");
    hData_rec->GetXaxis()->SetTitle(title+Reg);
    hData_subtracted->SetMarkerColor(4);
    hData_subtracted->SetLineStyle(2);
    hData_subtracted->SetLineColor(4);
    hData_subtracted->SetLineWidth(2);
    hData_subtracted->Draw("SAME");    

    TLegend *leg1 = new TLegend(x1,y1,x2,y2);
    leg1->SetBorderSize(1); 
    leg1->SetFillColor(29);
    leg1->AddEntry(hData_rec,"Data","pl");
    leg1->AddEntry(hData_subtracted,"Subtracted Data","pl");
    leg1->Draw();
    recVsSub_bef->SaveAs(output+varOutput+"recVsSub_bef_"+varToUnfold+"_"+varReg+".jpg");
    recVsSub_bef->SaveAs(output+varOutput+"recVsSub_bef_"+varToUnfold+"_"+varReg+".pdf");    


    //===== rec data Vs subtracted data after binning opt
    recVsSub_aft= new TCanvas("recVSbkg_aft_"+varToUnfold+"_"+varReg,"rec_Vs_Subtracted_data_"+varToUnfold+"_"+varReg,1000,1000);
    hnewData_rec->SetLineWidth(2);
    hnewData_rec->Draw();
    hnewData_rec->SetTitle("Rec Data Vs Subtracted Data");
    hnewData_rec->GetXaxis()->SetTitle(title+Reg);
    hnewData_subtracted->SetMarkerColor(4);
    hnewData_subtracted->SetLineColor(4);
    hnewData_subtracted->SetLineStyle(2);
    hnewData_subtracted->SetLineWidth(2);
    hnewData_subtracted->Draw("SAME");    

    TLegend *leg2 = new TLegend(x1,y1,x2,y2);
    leg2->SetBorderSize(1); 
    leg2->SetFillColor(29);
    leg2->AddEntry(hnewData_rec,"Data","pl");
    leg2->AddEntry(hnewData_subtracted,"Subtracted Data","pl");
    leg2->Draw();

    recVsSub_aft->SaveAs(output+varOutput+"recVsSub_aft_"+varToUnfold+"_"+varReg+".jpg");
    recVsSub_aft->SaveAs(output+varOutput+"recVsSub_aft_"+varToUnfold+"_"+varReg+".pdf");    
   
    //==== Unfolded data vs hGen
    unfoldedVsGen = new TCanvas("GenVsUnfold_"+varToUnfold+"_"+varReg, "GenVsUnfolded_"+varToUnfold+"_"+varReg,1000,1000);
    hData_unfolded->Scale(1./ hData_unfolded->Integral("width") );
    hnewMC_gen->Scale(1./ hnewMC_gen->Integral("width") );
    hData_unfolded->SetLineWidth(2);
    hData_unfolded->Draw();
    hData_unfolded->SetTitle("Unfolded data Vs Gen MC");
    hData_unfolded->GetXaxis()->SetTitle(title+Reg);
    hnewMC_gen->SetMarkerColor(4);
    hnewMC_gen->SetLineStyle(2);
    hnewMC_gen->SetLineColor(4);
    hnewMC_gen->SetLineWidth(2);
    hnewMC_gen->Draw("SAME");

    TLegend *leg3 = new TLegend(x1,y1,x2,y2);
    leg3->SetBorderSize(1);
    leg3->SetFillColor(29);
    leg3->AddEntry(hData_unfolded,"Unfolded Data","pl");
    leg3->AddEntry(hnewMC_gen,"MC Gen-Level","pl");
    leg3->Draw();

    unfoldedVsGen->SaveAs(output+varOutput+"UnfoldedVsGen_"+varToUnfold+"_"+varReg+".jpg");
    unfoldedVsGen->SaveAs(output+varOutput+"UnfoldedVsGen_"+varToUnfold+"_"+varReg+".pdf");    

    //==== Unfolded data subtracted vs hGen
    unfoldedSubVsGen = new TCanvas("GenVsUnfoldSub_"+varToUnfold+"_"+varReg, "GenVsUnfoldedSub_"+varToUnfold+"_"+varReg,1000,1000);
    hData_sub_unfolded->Scale(1./ hData_sub_unfolded->Integral("width") );
    hnewMC_gen->Scale(1./ hnewMC_gen->Integral("width") );
    hData_sub_unfolded->SetLineWidth(2);
    hData_sub_unfolded->Draw();
    hData_sub_unfolded->SetTitle("Unfolded Subtracted Data Vs Gen MC");
    hData_sub_unfolded->GetXaxis()->SetTitle(title+Reg);
    hnewMC_gen->SetMarkerColor(4);
    hnewMC_gen->SetLineColor(4);
    hnewMC_gen->SetLineStyle(2);
    hnewMC_gen->SetLineWidth(2);
    hnewMC_gen->Draw("SAME");
    
    TLegend *leg4 = new TLegend(x1,y1,x2,y2);
    leg4->SetBorderSize(1); 
    leg4->SetFillColor(29);
    leg4->AddEntry(hnewMC_gen,"MC Gen-Level","pl");
    leg4->AddEntry(hData_sub_unfolded,"Unfolded Subtracted Data","pl");
    leg4->Draw();

    unfoldedSubVsGen->SaveAs(output+varOutput+"UnfoldedSubVsGen_"+varToUnfold+"_"+varReg+".jpg");
    unfoldedSubVsGen->SaveAs(output+varOutput+"UnfoldedSubVsGen_"+varToUnfold+"_"+varReg+".pdf");    

    //==== Unfolded data vs hGen vs UnfoldedRec
    unfoldedVsRecVsGen = new TCanvas("GenVsSubUnfoldVsRecUnfold_"+varToUnfold+"_"+varReg, "GenVsSubUnfoldVsRecUnfold_"+varToUnfold+"_"+varReg,1000,1000);
    hData_sub_unfolded->Scale(1./ hData_unfolded->Integral("width") );
    hnewMC_gen->Scale(1./ hnewMC_gen->Integral("width") );
    hRec_unfolded->Scale(1./ hRec_unfolded->Integral("width") );
    hData_sub_unfolded->SetLineWidth(2);
    hData_sub_unfolded->Draw();
    hData_sub_unfolded->GetXaxis()->SetTitle(title+Reg);
    hData_sub_unfolded->SetTitle("Unfolded Sub data Vs Unfolded Rec MC Vs Gen MC");
    hnewMC_gen->SetLineWidth(2);
    hnewMC_gen->SetLineStyle(2);
    hnewMC_gen->SetMarkerColor(4);
    hnewMC_gen->SetLineColor(4);
    hnewMC_gen->Draw("SAME");
    hRec_unfolded->SetLineWidth(2);
    hRec_unfolded->SetMarkerColor(6);
    hRec_unfolded->SetLineStyle(3);
    hRec_unfolded->SetLineColor(6);
    hRec_unfolded->Draw("SAME");

    TLegend *leg5 = new TLegend(x1,y1,x2,y2);
    leg5->SetBorderSize(1); 
    leg5->SetFillColor(29);
    leg5->AddEntry(hData_sub_unfolded,"Unfolded Subtracted Data","pl");
    leg5->AddEntry(hnewMC_gen,"MC Gen-Level","pl");
    leg5->AddEntry(hRec_unfolded,"Unfolded MC Rec-Level","pl");
    leg5->Draw();

    unfoldedVsRecVsGen->SaveAs(output+varOutput+"UnfoldedDataVsRecMCVsGen_"+varToUnfold+"_"+varReg+".jpg");
    unfoldedVsRecVsGen->SaveAs(output+varOutput+"UnfoldedDataVsRecMCVsGen_"+varToUnfold+"_"+varReg+".pdf");    
  
    //==== Data sub Folded Vs Data Sub
    foldedVsSubData = new TCanvas("FoldedVsSubData_"+varToUnfold+"_"+varReg, "FoldedVsSubData_"+varToUnfold+"_"+varReg,1000,1000);
    hData_sub_folded->Scale(1./ hData_sub_folded->Integral("width"));
    hnewData_subtracted->Scale(1./ hnewData_subtracted->Integral("width"));
    hnewData_subtracted->SetMarkerColor(4);
    hnewData_subtracted->GetXaxis()->SetTitle(title+Reg);
    hnewData_subtracted->SetTitle("Folded Subtracted Data Vs Subtracted Data");
    hnewData_subtracted->SetLineColor(4);
    hnewData_subtracted->SetLineWidth(2);
    hnewData_subtracted->Draw();
    hData_sub_folded->SetLineStyle(2);
    hData_sub_folded->SetLineWidth(2);
    hData_sub_folded->Draw("HIST SAME");

    TLegend *leg6 = new TLegend(x1,y1,x2,y2);
    leg6->SetBorderSize(1);
    leg6->SetFillColor(29);
    leg6->AddEntry(hData_sub_folded,"Folded Subtracted Data","pl");
    leg6->AddEntry(hnewData_subtracted,"Subtracted Data","pl");
    leg6->Draw();

    foldedVsSubData->SaveAs(output+varOutput+"foldedVsSubData_"+varToUnfold+"_"+varReg+".jpg");
    foldedVsSubData->SaveAs(output+varOutput+"foldedVsSubData_"+varToUnfold+"_"+varReg+".pdf");
  }
}
