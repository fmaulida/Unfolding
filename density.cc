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
#include "TUnfoldSys.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TSVDUnfold.h"
#endif
#include "TUnfoldDensity.h"
#include "TChain.h"

float computeEfficiency(TH2 *h,int firstbin,int secondbin)
{
  float firstbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(firstbin) );
  float secondbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(secondbin) );
  float ngen_reco=h->Integral(firstbin,secondbin,firstbinReco,secondbinReco);
  float ngen=h->Integral();
  return ngen_reco/ngen;
}

float computePurity(TH2 *h,int firstbin,int secondbin)
{
  //  float pur(0);

  float firstbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(firstbin) );
  float secondbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(secondbin) );
						
  float ngen_reco=h->Integral(firstbin,secondbin,firstbinReco,secondbinReco);
  float ngen=h->Integral(firstbin,secondbin,1,h->GetYaxis()->GetNbins());

  return ngen>0 ? ngen_reco/ngen : 0;
}

float computeStability(TH2 *h,int firstbin,int secondbin)
{
  // float stab(0);

  float firstbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(firstbin) );
  float secondbinReco=h->GetYaxis()->FindBin( h->GetXaxis()->GetBinCenter(secondbin) );
						
  float ngen_reco=h->Integral(firstbin,secondbin,firstbinReco,secondbinReco);
  float nreco=h->Integral(1,h->GetXaxis()->GetNbins(),firstbinReco,secondbinReco);

  return nreco>0 ? ngen_reco/nreco : 0;
}

void density()
{
  //===========defined variable input==========================

  TString varToUnfold, varOutput, title;
  Int_t reg = 4,nbinsunfolded=100;
  TString varReg[4]={"away","tran","tow","inc"};
  Float_t xmin,xmax;
  Int_t TowRecBins=0, TowGenBins=0, TranRecBins=0, TranGenBins=0, AwayRecBins=0, AwayGenBins=0; 
  Float_t TowRebinRec[40]={ }; Float_t TowRebinGen[40]={ }, TranRebinRec[40]={ }, TranRebinGen[40]={ };
  Float_t AwayRebinRec[40]={ }, AwayRebinGen[40]={ };
  TString input = "/afs/cern.ch/work/f/fmaulida/CMSSW_7_2_0_pre3/src/sample/MC/";
  TString output = "/afs/cern.ch/work/f/fmaulida/13TeV/CMSSW_8_0_5_patch1/src/Unfolding/";
  Float_t rec_inc=0,gen_inc=0,rec_away=0,gen_away=0,rec_tran=0,gen_tran=0,rec_tow=0,gen_tow=0;
   
  cout<<"varToUnfold:"<<endl;
  cin>>varToUnfold;
  
  cout<<"Output:"<<endl;
  cin>>varOutput;

  if(varToUnfold=="avgptflux"){
    xmin=0;
    xmax=6;
    title = "<#sum p_{T} >";
  }
  if(varToUnfold=="ptflux"){
    xmin=0;
    xmax=100;
    title = "#sum p_{T}";
  }
  if(varToUnfold=="nch"){
    xmin=0;
    xmax=70;
    title = "N_{ch}";
  }


  //use N x 2N matrices
  int nbinsreconstructed=2*nbinsunfolded;

   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);

   TH1D *hGen[reg];
   TH1D *hRec[reg];
   TH2D *hGenRec[reg][reg];
   TH1D *hPurity[reg];
   TH1D *hStability[reg];
   TH2D *hnewGenRec[reg][reg];
   TH1D *hnewPurity[reg];
   TH1D *hnewStability[reg];
   TH1D *hnewRec[reg];
   TH1D *hnewGen[reg];
   TH1D *hnewDRec[reg];
   TH1D *hDRec[reg];
   TH1F *hTune_Zlep[reg];
   TH1F *hTuneZlep[reg];
   TH1F *hnewTune_Zlep[reg];
   TH1F *hTune_noCR[reg];
   TH1F *hTunenoCR[reg];
   TH1F *hnewTune_noCR[reg];
   TH1F *hTune_noMPI[reg];
   TH1F *hTunenoMPI[reg];
   TH1F *hnewTune_noMPI[reg];
   TH1F *hTune_UEP11[reg];
   TH1F *hTuneUEP11[reg];
   TH1F *hnewTune_UEP11[reg];

   //variables: generated, reconstructed, event weight                                                                 
   Float_t recNchAway, genNchAway, recNchTran, genNchTran, recNchTow, genNchTow;
   Float_t recPtfluxAway, genPtfluxAway, recPtfluxTran, genPtfluxTran, recPtfluxTow, genPtfluxTow; 
   Float_t recAvgPtfluxAway, genAvgPtfluxAway, recAvgPtfluxTran, genAvgPtfluxTran, recAvgPtfluxTow, genAvgPtfluxTow;
   Float_t weight;

   //readout signal expectations and fill response matrix                                                              
   TChain *signal  = new TChain("dataAnalyzer/ue");
   signal->Add(input+"MC8TeV_TT_Z2star_powheg_pythia_0.root");
   signal->Add(input+"MC8TeV_TT_Z2star_powheg_pythia_1.root");
   signal->Add(input+"MC8TeV_TT_Z2star_powheg_pythia_2.root");
   signal->Add(input+"MC8TeV_TT_Z2star_powheg_pythia_3.root");
   signal->Add(input+"MC8TeV_TT_Z2star_powheg_pythia_4.root");
   
   signal->SetBranchAddress("rec_nch_away", &recNchAway);
   signal->SetBranchAddress("rec_nch_tran", &recNchTran);
   signal->SetBranchAddress("rec_nch_tow", &recNchTow);
   
   signal->SetBranchAddress("rec_ptflux_away", &recPtfluxAway);
   signal->SetBranchAddress("rec_ptflux_tran", &recPtfluxTran);
   signal->SetBranchAddress("rec_ptflux_tow", &recPtfluxTow);

   signal->SetBranchAddress("rec_avgptflux_away", &recAvgPtfluxAway);
   signal->SetBranchAddress("rec_avgptflux_tran", &recAvgPtfluxTran);
   signal->SetBranchAddress("rec_avgptflux_tow", &recAvgPtfluxTow);

   signal->SetBranchAddress("gen_nch_away", &genNchAway);
   signal->SetBranchAddress("gen_nch_tran", &genNchTran);
   signal->SetBranchAddress("gen_nch_tow", &genNchTow);

   signal->SetBranchAddress("gen_ptflux_away", &genPtfluxAway);
   signal->SetBranchAddress("gen_ptflux_tran", &genPtfluxTran);
   signal->SetBranchAddress("gen_ptflux_tow", &genPtfluxTow);

   signal->SetBranchAddress("gen_avgptflux_away", &genAvgPtfluxAway);
   signal->SetBranchAddress("gen_avgptflux_tran", &genAvgPtfluxTran);
   signal->SetBranchAddress("gen_avgptflux_tow", &genAvgPtfluxTow);

   signal->SetBranchAddress("weight", &weight);

   for(Int_t iReg=0; iReg!=reg; ++iReg){
     hGen[iReg] = new TH1D("Gen_"+varToUnfold+"_"+varReg[iReg], "Gen_"+varToUnfold+"_"+varReg[iReg],nbinsunfolded, xmin,xmax);
     hGen[iReg]->Sumw2();
     hRec[iReg] = new TH1D("Rec_"+varToUnfold+"_"+varReg[iReg], "Rec_"+varToUnfold+"_"+varReg[iReg], nbinsunfolded, xmin, xmax);
     hRec[iReg]->Sumw2();
     for(Int_t jReg=0; jReg!=reg; ++jReg){
       hGenRec[iReg][jReg] = new TH2D("GenRec_"+varToUnfold+"_"+varReg[iReg]+"_"+varReg[jReg], "Migration_Matrix_"+varToUnfold+
				      "_"+varReg[iReg]+"_"+varReg[jReg]+";gen;rec",nbinsunfolded, xmin, xmax,nbinsreconstructed,xmin,xmax);
       hGenRec[iReg][jReg]->Sumw2();
     }
     
     //Histogram for Purity
     hPurity[iReg] = new TH1D("Purity_"+varToUnfold+"_"+varReg[iReg],"Purity_"+varToUnfold+"_"+varReg[iReg],nbinsunfolded,xmin,xmax);
     hPurity[iReg]->Sumw2();
     hStability[iReg] = new TH1D("Stability_"+varToUnfold+"_"+varReg[iReg],"Stability_"+varToUnfold+"_"+varReg[iReg],nbinsunfolded,xmin,xmax);
     hStability[iReg]->Sumw2();
   }
   
   // Load for Response Matrix                                                                                         
   for (Int_t i= 0; i<signal->GetEntries(); i++){
     signal->GetEntry(i);
     if(varToUnfold=="nch")
       {
	 rec_inc=recNchAway+recNchTow+recNchTran;
	 gen_inc=genNchAway+genNchTow+genNchTran;
	 rec_away=recNchAway;
	 gen_away=genNchAway;
	 
	 rec_tow=recNchTow;
	 gen_tow=genNchTow;
	 
	 rec_tran=recNchTran;
	 gen_tran=genNchTran;
       }
     if(varToUnfold=="ptflux")
       {
	 rec_inc=recPtfluxAway+recPtfluxTow+recPtfluxTran;
	 gen_inc=genPtfluxAway+genPtfluxTow+genPtfluxTran;
	 
	 rec_away=recPtfluxAway;
	 gen_away=genPtfluxAway;
	 
	 rec_tow=recPtfluxTow;
	 gen_tow=genPtfluxTow;
	 
	 rec_tran=recPtfluxTran;
	 gen_tran=genPtfluxTran;
       }
     if(varToUnfold=="avgptflux")
       {
	 //avgptfux = ptflux/nch
	 rec_inc=(recPtfluxAway+recPtfluxTow+recPtfluxTran)/(recNchAway+recNchTow+recNchTran);
	 gen_inc=(genPtfluxAway+genPtfluxTow+genPtfluxTran)/(genNchAway+genNchTow+genNchTran);
	 
	 rec_away=recPtfluxAway/recNchAway;
	 gen_away=genPtfluxAway/genNchAway;
	 
	 rec_tow=recPtfluxTow/recNchTow;
	 gen_tow=genPtfluxTow/genNchTow;
	 	 
	 rec_tran=recPtfluxTran/recNchTran;
	 gen_tran=genPtfluxTran/genNchTran;
       }
     
     for(Int_t iReg=0; iReg!=reg; ++iReg){
       if(iReg==0){
	 hRec[iReg]->Fill(rec_away,weight);                                                                                                                            
	 hGen[iReg]->Fill(gen_away,weight);                                                                                                                           
	 for(Int_t jReg=0; jReg!=reg; ++jReg){
	   if(jReg==0){
	   hGenRec[iReg][jReg]->Fill(gen_away,rec_away,weight);
	   }
	   if(jReg==1){
	     hGenRec[iReg][jReg]->Fill(gen_away,rec_tran,weight);
           }
	   if(jReg==2){
	     hGenRec[iReg][jReg]->Fill(gen_away,rec_tow,weight);
           }
	   if(jReg==3){
	     hGenRec[iReg][jReg]->Fill(gen_away,rec_inc,weight);
           }
	 }
       }
       if(iReg==1){
	 hRec[iReg]->Fill(rec_tran,weight);
	 hGen[iReg]->Fill(gen_tran,weight);
	 for(Int_t jReg=0; jReg!=reg; ++jReg){
           if(jReg==0){
	     hGenRec[iReg][jReg]->Fill(gen_tran,rec_away,weight);
           }
           if(jReg==1){
             hGenRec[iReg][jReg]->Fill(gen_tran,rec_tran,weight);
           }
           if(jReg==2){
             hGenRec[iReg][jReg]->Fill(gen_tran,rec_tow,weight);
           }
           if(jReg==3){
             hGenRec[iReg][jReg]->Fill(gen_tran,rec_inc,weight);
           }
         }
       }
       if(iReg==2){
	 hRec[iReg]->Fill(rec_tow,weight);
	 hGen[iReg]->Fill(gen_tow,weight);
	 for(Int_t jReg=0; jReg!=reg; ++jReg){
           if(jReg==0){
             hGenRec[iReg][jReg]->Fill(gen_tow,rec_away,weight);
           }
           if(jReg==1){
             hGenRec[iReg][jReg]->Fill(gen_tow,rec_tran,weight);
           }
           if(jReg==2){
	     hGenRec[iReg][jReg]->Fill(gen_tow,rec_tow,weight);
           }
           if(jReg==3){
             hGenRec[iReg][jReg]->Fill(gen_tow,rec_tow,weight);
           }
         }
       }
       if(iReg==3){
	 hRec[iReg]->Fill(rec_inc,weight);
	 hGen[iReg]->Fill(gen_inc,weight);
	 for(Int_t jReg=0; jReg!=reg; ++jReg){
           if(jReg==0){
             hGenRec[iReg][jReg]->Fill(gen_inc,rec_away,weight);
           }
           if(jReg==1){
             hGenRec[iReg][jReg]->Fill(gen_inc,rec_tran,weight);
           }
           if(jReg==2){
             hGenRec[iReg][jReg]->Fill(gen_inc,rec_tow,weight);
           }
	   if(jReg==3){
             hGenRec[iReg][jReg]->Fill(gen_inc,rec_inc,weight);
           }
         }
       }
     }
   } 
       
   
   //==================================================Purity and Stability====================================================   
   //=======================================calculate purity and stability  fix bin width=========================================
   for(int iRec = 1; iRec <=nbinsunfolded; iRec++){                                                                                                             
     for(Int_t iReg=0; iReg!=reg; ++iReg){                                                  

       float purity = computePurity(hGenRec[iReg][iReg],iRec,iRec);                                                                                                    
       float stability = computeStability(hGenRec[iReg][iReg],iRec,iRec); 
             
       hPurity[iReg]->SetBinContent(iRec,purity);                                                                                                                     
       hStability[iReg]->SetBinContent(iRec,stability);                                                                                
     }
   }                                        
   
   
   //==================================rebinning if Purity and Stability >50%================================================
   //loop over generator level bins
   float pur(0), stab(0), eff(0);
   // std::vector<float> newGenBins,newRecBins,newPur,newStab;
   for(Int_t iReg=0; iReg!=reg; ++iReg){//<==========Loop for The number of regions
     std::vector<float> newGenBins,newRecBins,newPur,newStab;
     const int nxBins = hRec[iReg]->GetNbinsX();
     newGenBins.push_back(hGen[iReg]->GetXaxis()->GetBinLowEdge(1));
     newRecBins.push_back(hRec[iReg]->GetXaxis()->GetBinLowEdge(1));
       
     for(int xbin=1; xbin<=nxBins;)
     {
	 int firstbin=xbin;
	 int secondbin=xbin-1;
	 
	 //	   cout << firstbin << endl;
	 do
	 {
	   secondbin++;
	     pur=computePurity(hGenRec[iReg][iReg],firstbin,secondbin);
	     stab=computeStability(hGenRec[iReg][iReg],firstbin,secondbin);
	     eff=computeEfficiency(hGenRec[iReg][iReg],firstbin,secondbin);
	     
	     //  cout << "\t" << secondbin << " " << pur << " " << stab << endl;     
	 }
	 while( (eff<0.01 || pur<0.5 || stab<0.5) && secondbin<=nxBins);//<====== Purity and Stability should be >50%, and efficiency >5%
	 //	 while( (pur<0.5 || stab<0.5) && secondbin<=nxBins);//<====== Purity and Stability should be >50%, and efficiency >5%
	 
	 //new bin found
	 float xdown=hGen[iReg]->GetXaxis()->GetBinLowEdge(firstbin);
	 float xup=hGen[iReg]->GetXaxis()->GetBinUpEdge(secondbin);
	 newGenBins.push_back(xup);
	 
	 newRecBins.push_back(0.5*(xup+xdown));
	 newRecBins.push_back(xup);
	 
	 //new Purity with new Binning
	 newPur.push_back(pur);
	 newStab.push_back(stab);
	 
	 //move to next bin
	 xbin=secondbin+1;
      
     }   

     //===================================Definning new bins for new Migration Matrix==============================================
     //rebinning the 2D hnewDetRes
     size_t nrecBins(newRecBins.size());
     Double_t rebinRec[nrecBins];
     for(size_t i=0; i<=nrecBins; i++){
       rebinRec[i]=newRecBins[i];
       if(iReg==0){
         AwayRecBins=newRecBins.size();
         AwayRebinRec[i]=newRecBins[i];
       }
       if(iReg==1){
         TranRecBins=newRecBins.size();
         TranRebinRec[i]=newRecBins[i];
       }
       if(iReg==2){
         TowRecBins=newRecBins.size();
	 TowRebinRec[i]=newRecBins[i];
       }
     }

     size_t ngenBins(newGenBins.size());
     Double_t rebinGen[ngenBins];
     for(size_t i=0; i<=ngenBins; i++){
       rebinGen[i]=newGenBins[i];
       if(iReg==0){
	 AwayGenBins=newGenBins.size();
	 AwayRebinGen[i]=newGenBins[i];
       }
       if(iReg==1){
         TranGenBins=newGenBins.size();
         TranRebinGen[i]=newGenBins[i];
       }
       if(iReg==2){
         TowGenBins=newGenBins.size();
         TowRebinGen[i]=newGenBins[i];
       }
     }

     //=========================new migration matrix after rebinning==========================================
     // Filling new Migration matrix
     for(Int_t jReg=0; jReg!=4; ++jReg){
       hnewGenRec[iReg][jReg] = new TH2D("newGenRec_"+varToUnfold+"_"+varReg[iReg]+"_"+varReg[jReg],"new_Migration_Matrix_"+varToUnfold+"_"
					 +varReg[iReg]+"_"+varReg[jReg]+";gen;rec",ngenBins-1,rebinGen,nrecBins-1,rebinRec);
       hnewGenRec[iReg][jReg]->Sumw2();
       TAxis *xaxis = hGenRec[iReg][jReg]->GetXaxis();
       TAxis *yaxis = hGenRec[iReg][jReg]->GetYaxis();
       for(int j=0; j<=hGenRec[iReg][jReg]->GetNbinsY(); j++){
	 for(int i=0; i<=hGenRec[iReg][jReg]->GetNbinsX();i++){
	   //	   hnewGenRec[iReg][jReg]->SetBinContent(i,j,hGenRec[iReg][jReg]->GetBinContent(i,j));
	   hnewGenRec[iReg][jReg]->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),hGenRec[iReg][jReg]->GetBinContent(i,j));
	 }
       }
     }
     //===============================New Purity and Stability after rebinning==============================================
     //Filling new purity
     hnewPurity[iReg] =new TH1D("newPurity_"+varToUnfold+"_"+varReg[iReg], "new_Purity_"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     hnewPurity[iReg]->Sumw2();
     for(int i=1;i<=hnewPurity[iReg]->GetNbinsX();i++){
       hnewPurity[iReg]->SetBinContent(i,newPur[i-1]);
     }
     
     //Filling new Stability
     hnewStability[iReg] =new TH1D("newStability_"+varToUnfold+"_"+varReg[iReg], "new_Stability_"+varToUnfold+"_"+varReg[iReg],ngenBins-1,rebinGen);
     hnewStability[iReg]->Sumw2();
     for(Int_t i=1;i<=hnewStability[iReg]->GetNbinsX();i++){
       hnewStability[iReg]->SetBinContent(i,newStab[i-1]);
     }


     //Rebin method
     hnewRec[iReg]= new TH1D("newRec_"+varToUnfold+"_"+varReg[iReg],"new_rec"+varToUnfold+"_"+varReg[iReg], nrecBins-1,rebinRec);
     hnewRec[iReg]->Sumw2();
     hnewRec[iReg] = (TH1D*)hRec[iReg]->Rebin(nrecBins-1,"newRec_"+varToUnfold+"_"+varReg[iReg],rebinRec);

     hnewGen[iReg]= new TH1D("newGen_"+varToUnfold+"_"+varReg[iReg],"new_gen"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     hnewGen[iReg]->Sumw2();
     hnewGen[iReg] = (TH1D*)hGen[iReg]->Rebin(ngenBins-1,"newGen_"+varToUnfold+"_"+varReg[iReg],rebinGen);


     /*          
     //histogram for rec and gen in all region                                               
     hnewRec[iReg]= new TH1D("newRec_"+varToUnfold+"_"+varReg[iReg],"new_rec"+varToUnfold+"_"+varReg[iReg], nrecBins-1,rebinRec);
     hnewRec[iReg]->Sumw2();
     TAxis *recAxis = hRec[iReg]->GetXaxis();
     for(Int_t i=0; i<=hRec[iReg]->GetNbinsX(); i++){
       //hnewRec[iReg]->SetBinContent(i,hRec[iReg]->GetBinContent(i));
       hnewRec[iReg]->Fill(recAxis->GetBinCenter(i),hRec[iReg]->GetBinContent(i));
     }
    
     hnewGen[iReg]= new TH1D("newGen_"+varToUnfold+"_"+varReg[iReg],"new_gen"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     hnewGen[iReg]->Sumw2();
     TAxis *genAxis = hGen[iReg]->GetXaxis();
     for(Int_t i=0; i<=hGen[iReg]->GetNbinsX(); i++){
       //hnewGen[iReg]->SetBinContent(i,hGen[iReg]->GetBinContent(i));
       hnewGen[iReg]->Fill(genAxis->GetBinCenter(i),hGen[iReg]->GetBinContent(i));
     }
     */
   
     //==============Data===========================================================
     hDRec[iReg] = new TH1D("DataRec_"+varToUnfold+"_"+varReg[iReg], "Data_Rec_"+varToUnfold+"_"+varReg[iReg],nbinsunfolded, xmin, xmax);
     hDRec[iReg]->Sumw2();

     Float_t DrecNchAway, DrecNchTran, DrecNchTow;
     Float_t DrecPtfluxAway, DrecPtfluxTran, DrecPtfluxTow;
     Float_t DrecAvgPtfluxAway, DrecAvgPtfluxTran, DrecAvgPtfluxTow;
     Float_t Dweight; 
     Float_t Drec_inc=0,Drec_away=0,Drec_tran=0,Drec_tow=0;

     TChain *ch_data = new TChain("dataAnalyzer/ue");
     TString directory = "/afs/cern.ch/work/f/fmaulida/CMSSW_7_2_0_pre3/src/sample/data/";
     ch_data->Add(directory+"Data8TeV_MuEG2012A_0.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012A_1.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012B_0.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012B_1.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012B_2.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012B_3.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012B_4.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_0.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_1.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_2.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_3.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_4.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012C_5.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_0.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_1.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_2.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_3.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_4.root");
     ch_data->Add(directory+"Data8TeV_MuEG2012D_5.root");

     
     ch_data->SetBranchAddress("rec_nch_away", &DrecNchAway);
     ch_data->SetBranchAddress("rec_nch_tran", &DrecNchTran);
     ch_data->SetBranchAddress("rec_nch_tow", &DrecNchTow);
     
     ch_data->SetBranchAddress("rec_ptflux_away", &DrecPtfluxAway);
     ch_data->SetBranchAddress("rec_ptflux_tran", &DrecPtfluxTran);
     ch_data->SetBranchAddress("rec_ptflux_tow", &DrecPtfluxTow);
     
     ch_data->SetBranchAddress("rec_avgptflux_away", &DrecAvgPtfluxAway);
     ch_data->SetBranchAddress("rec_avgptflux_tran", &DrecAvgPtfluxTran);
     ch_data->SetBranchAddress("rec_avgptflux_tow", &DrecAvgPtfluxTow);
     
     ch_data->SetBranchAddress("weight", &Dweight);

     for (Int_t i= 0; i<ch_data->GetEntries(); i++){
       ch_data->GetEntry(i);
       if(varToUnfold=="nch")
	 {
	   Drec_inc=DrecNchAway+DrecNchTow+DrecNchTran;
	   Drec_away=DrecNchAway;
	   Drec_tow=DrecNchTow;
	   Drec_tran=DrecNchTran;
	 }
       if(varToUnfold=="ptflux")
	 {
	   Drec_inc=DrecPtfluxAway+DrecPtfluxTow+DrecPtfluxTran;
	   Drec_away=DrecPtfluxAway;
	   Drec_tow=DrecPtfluxTow;
	   Drec_tran=DrecPtfluxTran;
	 }
       if(varToUnfold=="avgptflux")
	 {
	   //avgptfux = ptflux/nch                                                                                                                                     
	   Drec_inc=(DrecPtfluxAway+DrecPtfluxTow+DrecPtfluxTran)/(DrecNchAway+DrecNchTow+DrecNchTran);
	   Drec_away=DrecPtfluxAway/DrecNchAway;
	   Drec_tow=DrecPtfluxTow/DrecNchTow;
	   Drec_tran=DrecPtfluxTran/DrecNchTran;
	 }

       if(iReg==0){
	 hDRec[iReg]->Fill(Drec_away,Dweight);
       }
       if(iReg==1){
	 hDRec[iReg]->Fill(Drec_tran,Dweight);
       }
       if(iReg==2){
	 hDRec[iReg]->Fill(Drec_tow,Dweight);
       }
       if(iReg==3){
	 hDRec[iReg]->Fill(Drec_inc,Dweight);
       }
     }   

     hnewDRec[iReg] = new TH1D("newDataRec_"+varToUnfold+"_"+varReg[iReg], "new_Data_Rec_"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
     hnewDRec[iReg]->Sumw2();
     hnewDRec[iReg] = (TH1D*)hDRec[iReg]->Rebin(nrecBins-1,"newDRec_"+varToUnfold+"_"+varReg[iReg],rebinRec);

      //==================Tune================================
     TString input_1 = "/afs/cern.ch/work/f/fmaulida/13TeV/CMSSW_8_0_5_patch1/src/Unfolding/newLScan/";
     TFile *tune1 = new TFile(input_1+"out_1.root");
     TFile *tune2 = new TFile(input_1+"out_2.root");
     TFile *tune3 = new TFile(input_1+"out_3.root");
     TFile *tune4 = new TFile(input_1+"out_4.root");

     if(varToUnfold=="nch"){ 
       if(iReg==0){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chMult-r3"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]);
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chMult-r3"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chMult-r3"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chMult-r3"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==1){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chMult-r2"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]);
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chMult-r2"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chMult-r2"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chMult-r2"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==2){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chMult-r1");
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chMult-r1"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chMult-r1"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chMult-r1"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==3){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chMult-r0"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chMult-r0"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chMult-r0"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chMult-r0"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
     }
     
     if(varToUnfold=="ptflux"){
       if(iReg==0){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chSumPt-r3"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chSumPt-r3"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chSumPt-r3"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chSumPt-r3"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==1){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chSumPt-r2"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chSumPt-r2"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chSumPt-r2"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chSumPt-r2"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==2){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chSumPt-r1"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chSumPt-r1"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chSumPt-r1"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chSumPt-r1"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==3){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chSumPt-r0"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chSumPt-r0"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chSumPt-r0"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chSumPt-r0"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
     }
     if(varToUnfold=="avgptflux"){
       if(iReg==0){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chAvgPt-r3"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chAvgPt-r3"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chAvgPt-r3"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chAvgPt-r3"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==1){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chAvgPt-r2"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chAvgPt-r2"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chAvgPt-r2"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chAvgPt-r2"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==2){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chAvgPt-r1"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chAvgPt-r1"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chAvgPt-r1"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chAvgPt-r1"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
       if(iReg==3){
	 hTune_UEP11[iReg] = (TH1F*)tune1->Get("CMS_TOP_13_007/chAvgPt-r0"); 
	 hTune_UEP11[iReg]->SetName("UEP11_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noCR[iReg] = (TH1F*)tune2->Get("CMS_TOP_13_007/chAvgPt-r0"); 
	 hTune_noCR[iReg]->SetName("noCR_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_noMPI[iReg] = (TH1F*)tune3->Get("CMS_TOP_13_007/chAvgPt-r0"); 
	 hTune_noMPI[iReg]->SetName("noMPI_"+varToUnfold+"_"+varReg[iReg]); 
	 hTune_Zlep[iReg] = (TH1F*)tune4->Get("CMS_TOP_13_007/chAvgPt-r0"); 
	 hTune_Zlep[iReg]->SetName("Zlep_"+varToUnfold+"_"+varReg[iReg]); 
       }
     }

     hTuneZlep[iReg] = new TH1F("TuneZlep_"+varToUnfold+"_"+varReg[iReg], "TuneZlep_"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
     hTuneZlep[iReg] = (TH1F*)hTune_Zlep[iReg]->Rebin(nrecBins-1,"TuneZlep_"+varToUnfold+"_"+varReg[iReg],rebinRec);
     hTuneUEP11[iReg] = new TH1F("TuneUEP11_"+varToUnfold+"_"+varReg[iReg], "TuneUEP11_"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
     hTuneUEP11[iReg] = (TH1F*)hTune_UEP11[iReg]->Rebin(nrecBins-1,"TuneUEP11_"+varToUnfold+"_"+varReg[iReg],rebinRec);
     hTunenoCR[iReg] = new TH1F("TunenoCR_"+varToUnfold+"_"+varReg[iReg], "TunenoCR_"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
     hTunenoCR[iReg] = (TH1F*)hTune_noCR[iReg]->Rebin(nrecBins-1,"TunenoCR_"+varToUnfold+"_"+varReg[iReg],rebinRec);
     hTunenoMPI[iReg] = new TH1F("TunenoMPI_"+varToUnfold+"_"+varReg[iReg], "TunenoMPI_"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
     hTunenoMPI[iReg] = (TH1F*)hTune_noMPI[iReg]->Rebin(nrecBins-1,"TunenoMPI_"+varToUnfold+"_"+varReg[iReg],rebinRec);

     //================mergin Tune =============================
     
     hnewTune_UEP11[iReg]= new TH1F("tune_UEP11_"+varToUnfold+"_"+varReg[iReg],"tune_UEP11"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     //     hnewTune_Zlep[iReg]->Sumw2();
     hnewTune_UEP11[iReg] = (TH1F*)hTune_UEP11[iReg]->Rebin(ngenBins-1,"newUEP11_"+varToUnfold+"_"+varReg[iReg],rebinGen);

     hnewTune_noCR[iReg]= new TH1F("tune_noCR_"+varToUnfold+"_"+varReg[iReg],"tune_noCR"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     //     hnewTune_Zlep[iReg]->Sumw2();
     hnewTune_noCR[iReg] = (TH1F*)hTune_noCR[iReg]->Rebin(ngenBins-1,"newnoCR_"+varToUnfold+"_"+varReg[iReg],rebinGen);

     hnewTune_noMPI[iReg]= new TH1F("tune_noMPI_"+varToUnfold+"_"+varReg[iReg],"tune_noMPI"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     //     hnewTune_Zlep[iReg]->Sumw2();
     hnewTune_noMPI[iReg] = (TH1F*)hTune_noMPI[iReg]->Rebin(ngenBins-1,"newnoMPI_"+varToUnfold+"_"+varReg[iReg],rebinGen);

     hnewTune_Zlep[iReg]= new TH1F("tune_Zlep_"+varToUnfold+"_"+varReg[iReg],"tune_Zlep"+varToUnfold+"_"+varReg[iReg], ngenBins-1,rebinGen);
     //     hnewTune_Zlep[iReg]->Sumw2();
     hnewTune_Zlep[iReg] = (TH1F*)hTune_Zlep[iReg]->Rebin(ngenBins-1,"newZ2Lep_"+varToUnfold+"_"+varReg[iReg],rebinGen);
     

   }
   //=============================================Large Matrix==================================================================
   
   //====================================================3N axis for rec=========================================================
   Int_t nrec3Bins(TowRecBins+TranRecBins+AwayRecBins-2);//<=========================produce 27 binc for nch -2=25
   Float_t rebin3Rec[nrec3Bins];
   Int_t rec = 0;

   for(Int_t j=0; j<TowRecBins; ++j){
     rebin3Rec[rec++] = TowRebinRec[j];
   }
   for(Int_t j=1; j<TranRecBins; ++j){
     rebin3Rec[rec++] = TowRebinRec[TowRecBins-1]+TranRebinRec[j];
   }
   for(Int_t j=1; j<AwayRecBins; ++j){
     rebin3Rec[rec++] = TowRebinRec[TowRecBins-1]+TranRebinRec[TranRecBins-1]+AwayRebinRec[j];
   }
  
   //====Print xBins for 3 rec his
   for(Int_t bin=0; bin<nrec3Bins; ++bin){
     //cout<<"rebin3Rec\t"<<rebin3Rec[bin]<<endl;
   }

   //==============================================3N Gen-level================================================================
   Int_t ngen3Bins(TowGenBins+TranGenBins+AwayGenBins-2);//<=========================produce 27 binc for nch -2=25 
   //cout<<ngen3Bins<<endl;
   Float_t rebin3Gen[ngen3Bins];

   Int_t gen = 0;

   for(Int_t j=0; j<TowGenBins; ++j){
     rebin3Gen[gen++] = TowRebinGen[j];
   }
   for(Int_t j=1; j<TranGenBins; ++j){
     rebin3Gen[gen++] = TowRebinGen[TowGenBins-1]+TranRebinGen[j];
   }
   for(Int_t j=1; j<AwayGenBins; ++j){
     rebin3Gen[gen++] = TowRebinGen[TowGenBins-1]+TranRebinGen[TranGenBins-1]+AwayRebinGen[j];
   }

   //==========Print xBins for  Gen his
   for(Int_t bin=0; bin<ngen3Bins; ++bin){
     //cout<<"rebin3Gen\t"<<rebin3Gen[bin]<<endl;
   }

   
   //=============================================Define histogram for rec-level in three regions======================================================               
   TH1D *h3Rec = new TH1D("large_3rec_"+varToUnfold,"Global rec "+varToUnfold, nrec3Bins-1 ,rebin3Rec);
   h3Rec->Sumw2();
   for(Int_t i=1; i<TowRecBins; i++){
     h3Rec->SetBinContent(i,hnewRec[2]->GetBinContent(i));
   }
   for(Int_t i=1; i<TranRecBins; i++){
     h3Rec->SetBinContent(i+TowRecBins-1,hnewRec[1]->GetBinContent(i));
   }
   for(Int_t i=1; i<AwayRecBins; i++){
     h3Rec->SetBinContent(i+TowRecBins+TranRecBins-2,hnewRec[0]->GetBinContent(i));
   }

   TH1D *h3Gen = new TH1D("large_3gen_"+varToUnfold,"Global gen "+varToUnfold, ngen3Bins-1 ,rebin3Gen);
   h3Gen->Sumw2();   
   for(Int_t i=1; i<TowGenBins; i++){
     h3Gen->SetBinContent(i,hnewGen[2]->GetBinContent(i));
   }
   for(Int_t i=1; i<TranGenBins; i++){
     h3Gen->SetBinContent(i+TowGenBins-1,hnewGen[1]->GetBinContent(i));
   }
   for(Int_t i=1; i<AwayRecBins; i++){
     h3Gen->SetBinContent(i+TowGenBins+TranGenBins-2,hnewGen[0]->GetBinContent(i));
   }

   //================== Define data reconstructed for 3 regions=====================================================
   TH1D *h3DRec = new TH1D("large_3rec_data"+varToUnfold,"Global_rec_data_"+varToUnfold, nrec3Bins-1 ,rebin3Rec);
   h3DRec->Sumw2();
   for(Int_t i=1; i<TowRecBins; i++){
     h3DRec->SetBinContent(i,hnewDRec[2]->GetBinContent(i));
   }
   for(Int_t i=1; i<TranRecBins; i++){
     h3DRec->SetBinContent(i+TowRecBins-1,hnewDRec[1]->GetBinContent(i));
   }
   for(Int_t i=1; i<AwayRecBins; i++){
     h3DRec->SetBinContent(i+TowRecBins+TranRecBins-2,hnewDRec[0]->GetBinContent(i));
   }

   //==================Define global matrix==============================================================================
   TH2F *h3Mig = new TH2F("h3Mig", "Global Migration Matrix "+varToUnfold+";gen;rec", ngen3Bins-1, rebin3Gen, nrec3Bins-1, rebin3Rec);
   h3Mig->Sumw2();
   // TAxis *x3mig = h3Mig->GetXaxis();
   // TAxis *y3mig = h3Mig->GetYaxis();
  
   for(Int_t j=1; j<TowRecBins; j++){
     for(Int_t i=1; i<TowGenBins; i++){
       h3Mig->SetBinContent(i,j,hnewGenRec[2][2]->GetBinContent(i,j));
     }  
     for(Int_t i=1; i<TranGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins-1,j,hnewGenRec[1][2]->GetBinContent(i,j));
     }
     for(Int_t i=1; i<AwayGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins+TranGenBins-2,j,hnewGenRec[0][2]->GetBinContent(i,j));
     }
   }   
   for(Int_t j=1; j<TranRecBins; j++){
     for(Int_t i=1; i<TowGenBins; i++){
       h3Mig->SetBinContent(i,j+TowRecBins-1,hnewGenRec[2][1]->GetBinContent(i,j));
     }
     for(Int_t i=1; i<TranGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins-1,j+TowRecBins-1,hnewGenRec[1][1]->GetBinContent(i,j));
     }
     for(Int_t i=1; i<AwayGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins+TranGenBins-2,j+TowRecBins-1,hnewGenRec[0][1]->GetBinContent(i,j));
     }
   }
   for(Int_t j=1; j<AwayRecBins; j++){
     for(Int_t i=1; i<TowGenBins; i++){
       h3Mig->SetBinContent(i,j+TranRecBins+TowRecBins-2,hnewGenRec[2][0]->GetBinContent(i,j));
     }
     for(Int_t i=1; i<TranGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins-1,j+TranRecBins+TowRecBins-2,hnewGenRec[1][0]->GetBinContent(i,j));
     }
     for(Int_t i=1; i<AwayGenBins; i++){
       h3Mig->SetBinContent(i+TowGenBins+TranGenBins-2,j+TranRecBins+TowRecBins-2,hnewGenRec[0][0]->GetBinContent(i,j));
     }
   }
   
          
   //======================================== Do Unfolding starting from here =======================================
   //===========================define histogram for unfolded data===================================
   TH1 *data_unfolded[reg];
   TH1D *data_folded[reg];
   TGraph *rhoAvgGr[reg];
   TGraph *rhoMaxGr[reg];
   Float_t rhoAvg[reg];
   Float_t rhoMax[reg];
   Double_t tauUsed;
   Double_t scaleBias = 1;

   for(Int_t iReg=0; iReg!=reg; ++iReg){
     data_unfolded[iReg] = (TH1D*)hnewGen[iReg]->Clone("data_unfolded_"+varToUnfold+"_"+varReg[iReg]);
     // data_folded[iReg]=(TH1D*)hnewDRec[iReg]->Clone("data_folded");
     TUnfoldDensity unfold(hnewGenRec[iReg][iReg], TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeCurvature, TUnfold::kEConstraintArea);
     unfold.SetInput(hnewDRec[iReg],scaleBias);
     
     //============ scan tau regulation ================================================================================
     Double_t tauMin=3.0e-06;                             
     Double_t tauMax=2.0e-02;

     //Double_t tauMin=1.e-1;
     //Double_t tauMax=1.e+6;
     TSpline *logTauX,*logTauY;
     TSpline *rhoScan = 0;
     rhoAvgGr[iReg]=new TGraph();
     rhoMaxGr[iReg]=new TGraph();
     Int_t nScan=500;

     unfold.ScanTau(nScan,tauMin,tauMax,&rhoScan,TUnfoldDensity::kEScanTauRhoAvg);
     int iBest = unfold.ScanTau(nScan,tauMin,tauMax,&rhoScan,TUnfoldDensity::kEScanTauRhoAvg);

     // create graphs with one point to visualize best choice of tau
     Double_t t[1],rho[1],x[1],y[1];
     rhoScan->GetKnot(iBest,t[0],rho[0]);
     TGraph *bestRho=new TGraph(1,t,rho);
     Double_t *tAll=new Double_t[nScan],*rhoAll=new Double_t[nScan];
     for(Int_t i=0;i<nScan;i++) {
       rhoScan->GetKnot(i,tAll[i],rhoAll[i]);
       cout<<i<<"  "<<tAll[i]<<"  "<<rhoAll[i]<<endl;
     }

     cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  TAU,iBest --> "<<unfold.GetTau()<<"  "<<iBest<<endl;


     TGraph *knots=new TGraph(nScan,tAll,rhoAll);
     cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
	 <<" / "<<unfold.GetNdf()<<"\n";
   
     TCanvas *tcanv = new TCanvas();
     tcanv->cd();
     rhoScan->Draw();
     rhoScan->SetTitle(";log_{10}(#tau);average(#rho_{i})");
     rhoScan->SetLineColor(kRed);
     bestRho->Draw("*");
     tcanv->SaveAs("knots_"+varToUnfold+"_"+varReg[iReg]+".png");
     
     /* Double_t rhomin=10;
     
     for(Double_t tau=tauMin; tau<=tauMax; tau*=1.1)
       {
         unfold.DoUnfold(tau,hnewDRec[iReg],1.0);
         data_unfolded[iReg]->Reset("ICE");
         unfold.GetOutput(data_unfolded[iReg]);
         rhoAvg[iReg]=unfold.GetRhoAvg();
       //              cout<<"rhoAvgVStau"<<"\t"<<rhoAvg[iReg]<<"\t"<<tau<<endl 
         if(rhomin > rhoAvg[iReg]){
           rhomin = rhoAvg[iReg];
           tauUse[iReg] = tau;
	 }
       
         rhoAvgGr[iReg]->SetPoint(rhoAvgGr[iReg]->GetN(),tau,rhoAvg[iReg]);
         rhoMax[iReg]=unfold.GetRhoMax();
         rhoMaxGr[iReg]->SetPoint(rhoMaxGr[iReg]->GetN(),tau,rhoMax[iReg]);
	 // cout<<"rhoavgVstau"<<"\t"<<rhoAvg[iReg]<<"\t"<<tau<<endl;
       }
     cout<<"rhominVstauUse"<<"\t"<<rhomin<<"\t"<<tauUse[iReg]<<endl;
     cout<<tauUse[iReg]<<endl;
     */

     tauUsed = unfold.GetTau();
     unfold.DoUnfold(tauUsed,hnewDRec[iReg],scaleBias);
     data_unfolded[iReg]=unfold.GetOutput("unfolded");

     cout<<"TUnfold version: "<<TUnfold::GetTUnfoldVersion()<<"\n";
       
     TFile *out[reg];
     out[iReg]= new TFile(output+varOutput+"out_"+varToUnfold+"_"+varReg[iReg]+".root", "RECREATE");
     
     hGen[iReg]->Write();
     hRec[iReg]->Write();
     hGenRec[iReg][iReg]->Write();
     hPurity[iReg]->Write();
     hStability[iReg]->Write();
     hnewGenRec[iReg][iReg]->Write();
     hnewPurity[iReg]->Write();
     hnewStability[iReg]->Write();
     hnewRec[iReg]->Write();
     hnewGen[iReg]->Write();
     hnewDRec[iReg]->Write();
     hDRec[iReg]->Write();
     hTune_Zlep[iReg]->Write();
     hTuneZlep[iReg]->Write();
     hTuneUEP11[iReg]->Write();
     hTunenoMPI[iReg]->Write();
     hTunenoCR[iReg]->Write();
     hnewTune_Zlep[iReg]->Write();
     hTune_noCR[iReg]->Write();
     hnewTune_noCR[iReg]->Write();
     hTune_noMPI[iReg]->Write();
     hnewTune_noMPI[iReg]->Write();
     hTune_UEP11[iReg]->Write();
     hnewTune_UEP11[iReg]->Write();
     //rhoMaxGr[iReg]->Write();
     //  rhoAvgGr[iReg]->Write();
     data_unfolded[iReg]->Write();
     out[iReg]->Close();
    
   }
}
