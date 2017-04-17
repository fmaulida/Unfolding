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

void Opt_binning()
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
  TString output = "/afs/cern.ch/work/f/fmaulida/13TeV/CMSSW_8_0_5_patch1/src/Unfolding/newlooks/";
  Float_t rec_inc=0,gen_inc=0,rec_away=0,gen_away=0,rec_tran=0,gen_tran=0,rec_tow=0,gen_tow=0;
  
  //use N x 2N matrices
  int nbinsreconstructed=2*nbinsunfolded;

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

  //===MC==========
  TH1D *hRec[reg];
  TH1D *hGen[reg];
  TH1D *hnewRec[reg];
  TH1D *hnewGen[reg];
  TH2D *hGenRec[reg][reg];
  TH1D *hPurity[reg];
  TH1D *hStability[reg];
  TH2D *hnewGenRec[reg][reg];
  TH1D *hnewPurity[reg];
  TH1D *hnewStability[reg];
  
  //===Data=
  TH1D *hDRec[reg];
  TH1D *hnewDRec[reg];
  TH1D *hData_subtracted[reg];
  TH1D *hnewData_subtracted[reg];
  
  //===Bkg====
  TH1D *hBkg_WW[reg];
  TH1D *hBkg_WZ[reg];
  TH1D *hBkg_ZZ[reg];
  TH1D *hBkg_singleT_tW[reg];
  TH1D *hBkg_singleTbar_tW[reg];
  
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
    hRec[iReg] = new TH1D("Rec_"+varToUnfold+"_"+varReg[iReg], "Rec_"+varToUnfold+"_"+varReg[iReg],nbinsunfolded, xmin, xmax);
    hRec[iReg]->Sumw2();
    for(Int_t jReg=0; jReg!=reg; ++jReg){
      hGenRec[iReg][jReg] = new TH2D("GenRec_"+varToUnfold+"_"+varReg[iReg]+"_"+varReg[jReg], "Migration_Matrix_"+varToUnfold+
				     "_"+varReg[iReg]+"_"+varReg[jReg]+";gen;rec",nbinsunfolded, xmin, xmax,nbinsreconstructed,xmin,xmax);
      hGenRec[iReg][jReg]->Sumw2();
    }
     
    //Histogram for Purity & Stability
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
  for(int iGen = 1; iGen <=nbinsunfolded; iGen++){ 
    for(Int_t iReg=0; iReg!=reg; ++iReg){                                                  

      float purity = computePurity(hGenRec[iReg][iReg],iGen,iGen);                                                
      float stability = computeStability(hGenRec[iReg][iReg],iGen,iGen);

      hPurity[iReg]->SetBinContent(iGen,purity);      
      hStability[iReg]->SetBinContent(iGen,stability);             
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
	 
	//   cout << firstbin << endl;
	 do
	   {
	     secondbin++;
	     pur=computePurity(hGenRec[iReg][iReg],firstbin,secondbin);
	     stab=computeStability(hGenRec[iReg][iReg],firstbin,secondbin);
	     eff=computeEfficiency(hGenRec[iReg][iReg],firstbin,secondbin);
	          
	     //  cout << "\t" << secondbin << " " << pur << " " << stab << endl;     
	   }
	 while( (eff<0.01 || pur<0.5 || stab<0.5) && secondbin<=nxBins);//<====== Purity and Stability should be >50%, and efficiency >5%
	 // while( (pur<0.5 || stab<0.5) && secondbin<=nxBins);//<====== Purity and Stability should be >50%, and efficiency >5%
	  
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
	  //   hnewGenRec[iReg][jReg]->SetBinContent(i,j,hGenRec[iReg][jReg]->GetBinContent(i,j));
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


    //==============Data Binning===========================================
    Float_t DrecNchAway, DrecNchTran, DrecNchTow;
    Float_t DrecPtfluxAway, DrecPtfluxTran, DrecPtfluxTow;
    Float_t DrecAvgPtfluxAway, DrecAvgPtfluxTran, DrecAvgPtfluxTow;
    Float_t Dweight;
    Float_t Drec_inc=0,Drec_away=0,Drec_tran=0,Drec_tow=0;

    TChain *ch_data = new TChain("dataAnalyzer/ue");
    TString directory = "/afs/cern.ch/work/f/fmaulida/CMSSW_7_2_0_pre3/src/sample/data/";              
    //TString directory = "/home/ami/pythia_root/data/";
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

    hDRec[iReg] = new TH1D("Rec_data_"+varToUnfold+"_"+varReg[iReg], "Rec_data_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hDRec[iReg]->Sumw2();

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
	  //avgptfux = ptflux/nch                                                                           \
                                                                                                              
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

    //==========Rebin Data rec================
    hnewDRec[iReg] = new TH1D("newRec_data"+varToUnfold+"_"+varReg[iReg], "new_Rec_data"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
    hnewDRec[iReg]->Sumw2();
    hnewDRec[iReg] = (TH1D*)hDRec[iReg]->Rebin(nrecBins-1,"new_Rec_data"+varToUnfold+"_"+varReg[iReg],rebinRec);
    
    
    //=================Bkg==========================================================================
    //===============WW===============================================================
    //----------------------------------------------Background _WW----------------------------------
    Float_t BrecNchAway, BrecNchTran, BrecNchTow;
    Float_t BrecPtfluxAway, BrecPtfluxTran, BrecPtfluxTow;
    Float_t BrecAvgPtfluxAway, BrecAvgPtfluxTran, BrecAvgPtfluxTow;
    Float_t Bweight; 
    Float_t Brec_inc=0,Brec_away=0,Brec_tran=0,Brec_tow=0;
    float lumi = 19701.0;

    TChain *bkg_WW = new TChain("dataAnalyzer/ue");
    bkg_WW->Add(input+"MC8TeV_WW_0.root");

    TFile* WW=new TFile(input+"MC8TeV_WW_0.root", "update");
    TVectorD *constValsWW = (TVectorD *)WW->Get("constVals");
    float norigEventsWW = (*constValsWW)[0];
    float xsecWW = (*constValsWW)[1];
    float xsecWeightWW = lumi*xsecWW/norigEventsWW;

    bkg_WW->SetBranchAddress("rec_nch_away", &BrecNchAway);
    bkg_WW->SetBranchAddress("rec_nch_tran", &BrecNchTran);
    bkg_WW->SetBranchAddress("rec_nch_tow", &BrecNchTow);

    bkg_WW->SetBranchAddress("rec_ptflux_away", &BrecPtfluxAway);
    bkg_WW->SetBranchAddress("rec_ptflux_tran", &BrecPtfluxTran);
    bkg_WW->SetBranchAddress("rec_ptflux_tow", &BrecPtfluxTow);

    bkg_WW->SetBranchAddress("rec_avgptflux_away", &BrecAvgPtfluxAway);
    bkg_WW->SetBranchAddress("rec_avgptflux_tran", &BrecAvgPtfluxTran);
    bkg_WW->SetBranchAddress("rec_avgptflux_tow", &BrecAvgPtfluxTow);

    bkg_WW->SetBranchAddress("weight", &Bweight);

    //histogram for Background                                                                                                                                                                                
    hBkg_WW[iReg]  = new TH1D("Bkg_WW_"+varToUnfold+"_"+varReg[iReg], "bkg_WW_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hBkg_WW[iReg]->Sumw2();
    for (Int_t i= 0; i<bkg_WW->GetEntries(); i++) {
      bkg_WW->GetEntry(i);
      float eventWeight(xsecWeightWW*Bweight);
      if(varToUnfold=="nch")
	{
	  Brec_inc=BrecNchAway+BrecNchTow+BrecNchTran;
	  Brec_away=BrecNchAway;
	  Brec_tow=BrecNchTow;
	  Brec_tran=BrecNchTran;
	}
      if(varToUnfold=="ptflux")
	{
	  Brec_inc=BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran;
	  Brec_away=BrecPtfluxAway;
	  Brec_tow=BrecPtfluxTow;
	  Brec_tran=BrecPtfluxTran;
	}
      if(varToUnfold=="avgptflux")
	{
	  //avgptfux = ptflux/nch                                                                                    
	  Brec_inc=(BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran)/(BrecNchAway+BrecNchTow+BrecNchTran);
	  Brec_away=BrecPtfluxAway/BrecNchAway;
	  Brec_tow=BrecPtfluxTow/BrecNchTow;
	  Brec_tran=BrecPtfluxTran/BrecNchTran;
	}
      if(iReg==0){
	hBkg_WW[iReg]->Fill(Brec_away,eventWeight);
      }
      if(iReg==1){
	hBkg_WW[iReg]->Fill(Brec_tran,eventWeight);
      }
      if(iReg==2){
	hBkg_WW[iReg]->Fill(Brec_tow,eventWeight);
      }
      if(iReg==3){
	hBkg_WW[iReg]->Fill(Brec_inc,eventWeight);
      }
    }

    //----------------------------------------------Background _WZ----------------------------------
    TChain *bkg_WZ = new TChain("dataAnalyzer/ue");
    bkg_WZ->Add(input+"MC8TeV_WZ_0.root");

    TFile* WZ=new TFile(input+"MC8TeV_WZ_0.root", "update");
    TVectorD *constValsWZ = (TVectorD *)WZ->Get("constVals");
    float norigEventsWZ = (*constValsWZ)[0];
    float xsecWZ = (*constValsWZ)[1];
    float xsecWeightWZ = lumi*xsecWZ/norigEventsWZ;

    bkg_WZ->SetBranchAddress("rec_nch_away", &BrecNchAway);
    bkg_WZ->SetBranchAddress("rec_nch_tran", &BrecNchTran);
    bkg_WZ->SetBranchAddress("rec_nch_tow", &BrecNchTow);

    bkg_WZ->SetBranchAddress("rec_ptflux_away", &BrecPtfluxAway);
    bkg_WZ->SetBranchAddress("rec_ptflux_tran", &BrecPtfluxTran);
    bkg_WZ->SetBranchAddress("rec_ptflux_tow", &BrecPtfluxTow);

    bkg_WZ->SetBranchAddress("rec_avgptflux_away", &BrecAvgPtfluxAway);
    bkg_WZ->SetBranchAddress("rec_avgptflux_tran", &BrecAvgPtfluxTran);
    bkg_WZ->SetBranchAddress("rec_avgptflux_tow", &BrecAvgPtfluxTow);

    bkg_WZ->SetBranchAddress("weight", &Bweight);

    //histogram for Background                                                                                                                                                                                
    hBkg_WZ[iReg]  = new TH1D("Bkg_WZ_"+varToUnfold+"_"+varReg[iReg], "bkg_WZ_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hBkg_WZ[iReg]->Sumw2();
    for (Int_t i= 0; i<bkg_WZ->GetEntries(); i++) {
      bkg_WZ->GetEntry(i);
      float eventWeight(xsecWeightWZ*Bweight);
      if(varToUnfold=="nch")
	{
	  Brec_inc=BrecNchAway+BrecNchTow+BrecNchTran;
	  Brec_away=BrecNchAway;
	  Brec_tow=BrecNchTow;
	  Brec_tran=BrecNchTran;
	}
      if(varToUnfold=="ptflux")
	{
	  Brec_inc=BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran;
	  Brec_away=BrecPtfluxAway;
	  Brec_tow=BrecPtfluxTow;
	  Brec_tran=BrecPtfluxTran;
	}
      if(varToUnfold=="avgptflux")
	{
	  //avgptfux = ptflux/nch                                                                                    
	  Brec_inc=(BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran)/(BrecNchAway+BrecNchTow+BrecNchTran);
	  Brec_away=BrecPtfluxAway/BrecNchAway;
	  Brec_tow=BrecPtfluxTow/BrecNchTow;
	  Brec_tran=BrecPtfluxTran/BrecNchTran;
	}
      if(iReg==0){
	hBkg_WZ[iReg]->Fill(Brec_away,eventWeight);
      }
      if(iReg==1){
	hBkg_WZ[iReg]->Fill(Brec_tran,eventWeight);
      }
      if(iReg==2){
	hBkg_WZ[iReg]->Fill(Brec_tow,eventWeight);
      }
      if(iReg==3){
	hBkg_WZ[iReg]->Fill(Brec_inc,eventWeight);
      }
    }

    //----------------------------------------------Background _ZZ----------------------------------
    TChain *bkg_ZZ = new TChain("dataAnalyzer/ue");
    bkg_ZZ->Add(input+"MC8TeV_ZZ_0.root");

    TFile* ZZ=new TFile(input+"MC8TeV_ZZ_0.root", "update");
    TVectorD *constValsZZ = (TVectorD *)ZZ->Get("constVals");
    float norigEventsZZ = (*constValsZZ)[0];
    float xsecZZ = (*constValsZZ)[1];
    float xsecWeightZZ = lumi*xsecZZ/norigEventsZZ;


    bkg_ZZ->SetBranchAddress("rec_nch_away", &BrecNchAway);
    bkg_ZZ->SetBranchAddress("rec_nch_tran", &BrecNchTran);
    bkg_ZZ->SetBranchAddress("rec_nch_tow", &BrecNchTow);

    bkg_ZZ->SetBranchAddress("rec_ptflux_away", &BrecPtfluxAway);
    bkg_ZZ->SetBranchAddress("rec_ptflux_tran", &BrecPtfluxTran);
    bkg_ZZ->SetBranchAddress("rec_ptflux_tow", &BrecPtfluxTow);

    bkg_ZZ->SetBranchAddress("rec_avgptflux_away", &BrecAvgPtfluxAway);
    bkg_ZZ->SetBranchAddress("rec_avgptflux_tran", &BrecAvgPtfluxTran);
    bkg_ZZ->SetBranchAddress("rec_avgptflux_tow", &BrecAvgPtfluxTow);

    bkg_ZZ->SetBranchAddress("weight", &Bweight);

    //histogram for Background                                                                                                                                                                                
    hBkg_ZZ[iReg]  = new TH1D("Bkg_ZZ_"+varToUnfold+"_"+varReg[iReg], "bkg_ZZ_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hBkg_ZZ[iReg]->Sumw2();
    for (Int_t i= 0; i<bkg_ZZ->GetEntries(); i++) {
      bkg_ZZ->GetEntry(i);
      float eventWeight(xsecWeightZZ*Bweight);
      if(varToUnfold=="nch")
	{
	  Brec_inc=BrecNchAway+BrecNchTow+BrecNchTran;
	  Brec_away=BrecNchAway;
	  Brec_tow=BrecNchTow;
	  Brec_tran=BrecNchTran;
	}
      if(varToUnfold=="ptflux")
	{
	  Brec_inc=BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran;
	  Brec_away=BrecPtfluxAway;
	  Brec_tow=BrecPtfluxTow;
	  Brec_tran=BrecPtfluxTran;
	}
      if(varToUnfold=="avgptflux")
	{
	  //avgptfux = ptflux/nch                                                                                    
	  Brec_inc=(BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran)/(BrecNchAway+BrecNchTow+BrecNchTran);
	  Brec_away=BrecPtfluxAway/BrecNchAway;
	  Brec_tow=BrecPtfluxTow/BrecNchTow;
	  Brec_tran=BrecPtfluxTran/BrecNchTran;
	}
      if(iReg==0){
	hBkg_ZZ[iReg]->Fill(Brec_away,eventWeight);
      }
      if(iReg==1){
	hBkg_ZZ[iReg]->Fill(Brec_tran,eventWeight);
      }
      if(iReg==2){
	hBkg_ZZ[iReg]->Fill(Brec_tow,eventWeight);
      }
      if(iReg==3){
	hBkg_ZZ[iReg]->Fill(Brec_inc,eventWeight);
      }
    }
         
    //----------------------------------------------Background _singleT_tW----------------------------------
    TChain *bkg_singleT_tW = new TChain("dataAnalyzer/ue");
    bkg_singleT_tW->Add(input+"MC8TeV_SingleT_tW.root");

    TFile* singleT_tW=new TFile(input+"MC8TeV_SingleT_tW.root", "update");
    TVectorD *constValssingleT_tW = (TVectorD *)singleT_tW->Get("constVals");
    float norigEventssingleT_tW = (*constValssingleT_tW)[0];
    float xsecsingleT_tW = (*constValssingleT_tW)[1];
    float xsecWeightsingleT_tW = lumi*xsecsingleT_tW/norigEventssingleT_tW;
     
    bkg_singleT_tW->SetBranchAddress("rec_nch_away", &BrecNchAway);
    bkg_singleT_tW->SetBranchAddress("rec_nch_tran", &BrecNchTran);
    bkg_singleT_tW->SetBranchAddress("rec_nch_tow", &BrecNchTow);

    bkg_singleT_tW->SetBranchAddress("rec_ptflux_away", &BrecPtfluxAway);
    bkg_singleT_tW->SetBranchAddress("rec_ptflux_tran", &BrecPtfluxTran);
    bkg_singleT_tW->SetBranchAddress("rec_ptflux_tow", &BrecPtfluxTow);

    bkg_singleT_tW->SetBranchAddress("rec_avgptflux_away", &BrecAvgPtfluxAway);
    bkg_singleT_tW->SetBranchAddress("rec_avgptflux_tran", &BrecAvgPtfluxTran);
    bkg_singleT_tW->SetBranchAddress("rec_avgptflux_tow", &BrecAvgPtfluxTow);

    bkg_singleT_tW->SetBranchAddress("weight", &Bweight);

    //histogram for Background                                                                                                                                                                                
    hBkg_singleT_tW[iReg]  = new TH1D("Bkg_singleT_tW_"+varToUnfold+"_"+varReg[iReg], "bkg_singleT_tW_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hBkg_singleT_tW[iReg]->Sumw2();
    for (Int_t i= 0; i<bkg_singleT_tW->GetEntries(); i++) {
      bkg_singleT_tW->GetEntry(i);
      float eventWeight(xsecWeightsingleT_tW*Bweight);
      if(varToUnfold=="nch")
	{
	  Brec_inc=BrecNchAway+BrecNchTow+BrecNchTran;
	  Brec_away=BrecNchAway;
	  Brec_tow=BrecNchTow;
	  Brec_tran=BrecNchTran;
	}
      if(varToUnfold=="ptflux")
	{
	  Brec_inc=BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran;
	  Brec_away=BrecPtfluxAway;
	  Brec_tow=BrecPtfluxTow;
	  Brec_tran=BrecPtfluxTran;
	}
      if(varToUnfold=="avgptflux")
	{
	  //avgptfux = ptflux/nch                                                                                    
	  Brec_inc=(BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran)/(BrecNchAway+BrecNchTow+BrecNchTran);
	  Brec_away=BrecPtfluxAway/BrecNchAway;
	  Brec_tow=BrecPtfluxTow/BrecNchTow;
	  Brec_tran=BrecPtfluxTran/BrecNchTran;
	}
      if(iReg==0){
	hBkg_singleT_tW[iReg]->Fill(Brec_away,eventWeight);
      }
      if(iReg==1){
	hBkg_singleT_tW[iReg]->Fill(Brec_tran,eventWeight);
      }
      if(iReg==2){
	hBkg_singleT_tW[iReg]->Fill(Brec_tow,eventWeight);
      }
      if(iReg==3){
	hBkg_singleT_tW[iReg]->Fill(Brec_inc,eventWeight);
      }
    }
     
    //----------------------------------------------Background _singleTbar_tW----------------------------------
    TChain *bkg_singleTbar_tW = new TChain("dataAnalyzer/ue");
    bkg_singleTbar_tW->Add(input+"MC8TeV_SingleTbar_tW.root");

    TFile* singleTbar_tW=new TFile(input+"MC8TeV_SingleTbar_tW.root", "update");
    TVectorD *constValssingleTbar_tW = (TVectorD *)singleTbar_tW->Get("constVals");
    float norigEventssingleTbar_tW = (*constValssingleTbar_tW)[0];
    float xsecsingleTbar_tW = (*constValssingleTbar_tW)[1];
    float xsecWeightsingleTbar_tW = lumi*xsecsingleTbar_tW/norigEventssingleTbar_tW;

    bkg_singleTbar_tW->SetBranchAddress("rec_nch_away", &BrecNchAway);
    bkg_singleTbar_tW->SetBranchAddress("rec_nch_tran", &BrecNchTran);
    bkg_singleTbar_tW->SetBranchAddress("rec_nch_tow", &BrecNchTow);

    bkg_singleTbar_tW->SetBranchAddress("rec_ptflux_away", &BrecPtfluxAway);
    bkg_singleTbar_tW->SetBranchAddress("rec_ptflux_tran", &BrecPtfluxTran);
    bkg_singleTbar_tW->SetBranchAddress("rec_ptflux_tow", &BrecPtfluxTow);

    bkg_singleTbar_tW->SetBranchAddress("rec_avgptflux_away", &BrecAvgPtfluxAway);
    bkg_singleTbar_tW->SetBranchAddress("rec_avgptflux_tran", &BrecAvgPtfluxTran);
    bkg_singleTbar_tW->SetBranchAddress("rec_avgptflux_tow", &BrecAvgPtfluxTow);

    bkg_singleTbar_tW->SetBranchAddress("weight", &Bweight);

    //histogram for Background                                                                                                                                                                                
    hBkg_singleTbar_tW[iReg]  = new TH1D("Bkg_singleTbar_tW_"+varToUnfold+"_"+varReg[iReg], "bkg_singleTbar_tW_"+varToUnfold+"_"+varReg[iReg], nbinsreconstructed, xmin, xmax);
    hBkg_singleTbar_tW[iReg]->Sumw2();
    for (Int_t i= 0; i<bkg_singleTbar_tW->GetEntries(); i++) {
      bkg_singleTbar_tW->GetEntry(i);
      float eventWeight(xsecWeightsingleTbar_tW*Bweight);
      if(varToUnfold=="nch")
	{
	  Brec_inc=BrecNchAway+BrecNchTow+BrecNchTran;
	  Brec_away=BrecNchAway;
	  Brec_tow=BrecNchTow;
	  Brec_tran=BrecNchTran;
	}
      if(varToUnfold=="ptflux")
	{
	  Brec_inc=BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran;
	  Brec_away=BrecPtfluxAway;
	  Brec_tow=BrecPtfluxTow;
	  Brec_tran=BrecPtfluxTran;
	}
      if(varToUnfold=="avgptflux")
	{
	  //avgptfux = ptflux/nch                                                                                    
	  Brec_inc=(BrecPtfluxAway+BrecPtfluxTow+BrecPtfluxTran)/(BrecNchAway+BrecNchTow+BrecNchTran);
	  Brec_away=BrecPtfluxAway/BrecNchAway;
	  Brec_tow=BrecPtfluxTow/BrecNchTow;
	  Brec_tran=BrecPtfluxTran/BrecNchTran;
	}
      if(iReg==0){
	hBkg_singleTbar_tW[iReg]->Fill(Brec_away,eventWeight);
      }
      if(iReg==1){
	hBkg_singleTbar_tW[iReg]->Fill(Brec_tran,eventWeight);
      }
      if(iReg==2){
	hBkg_singleTbar_tW[iReg]->Fill(Brec_tow,eventWeight);
      }
      if(iReg==3){
	hBkg_singleTbar_tW[iReg]->Fill(Brec_inc,eventWeight);
      }
    }
     
    //=====Data_subtracted==
    hData_subtracted[iReg]=(TH1D*)hDRec[iReg]->Clone("Data_sub_"+varToUnfold+"_"+varReg[iReg]);
    hData_subtracted[iReg]->Add(hBkg_WW[iReg],-1);
    hData_subtracted[iReg]->Add(hBkg_WZ[iReg],-1);
    hData_subtracted[iReg]->Add(hBkg_ZZ[iReg],-1);
    hData_subtracted[iReg]->Add(hBkg_singleT_tW[iReg],-1);
    hData_subtracted[iReg]->Add(hBkg_singleTbar_tW[iReg],-1);
    
    //==Rebin Data_subtracted===
    hnewData_subtracted[iReg]= new TH1D("newData_sub_"+varToUnfold+"_"+varReg[iReg], "new_data_sub"+varToUnfold+"_"+varReg[iReg],nrecBins-1,rebinRec);
    hnewData_subtracted[iReg]->Sumw2();
    hnewData_subtracted[iReg] = (TH1D*)hData_subtracted[iReg]->Rebin(nrecBins-1,"newData_sub_"+varToUnfold+"_"+varReg[iReg],rebinRec);
    
    //===Save MC========
    TFile *out_MC[reg];
    out_MC[iReg]= new TFile(output+varOutput+"out_MC_binning.root", "UPDATE");
     
    hGen[iReg]->Write();
    hnewGen[iReg]->Write();
    
    hRec[iReg]->Write();
    hnewRec[iReg]->Write();

    hGenRec[iReg][iReg]->Write();
    hnewGenRec[iReg][iReg]->Write();
    
    hPurity[iReg]->Write();
    hnewPurity[iReg]->Write();

    hStability[iReg]->Write();
    hnewStability[iReg]->Write();     
    
    out_MC[iReg]->Close();

    //=============Save Data======================
    TFile *out_Data[reg];
    out_Data[iReg]= new TFile(output+varOutput+"out_data_binning.root", "UPDATE");

    hDRec[iReg]->Write();
    hnewDRec[iReg]->Write();
    hData_subtracted[iReg]->Write();
    hnewData_subtracted[iReg]->Write();

    out_Data[iReg]->Close();


    //=====Save Bkg============
    TFile *out_Bkg[reg];
    out_Bkg[iReg]= new TFile(output+varOutput+"out_bkg_binning.root", "UPDATE");
    
    hBkg_WW[iReg]->Write();
    hBkg_WZ[iReg]->Write();
    hBkg_ZZ[iReg]->Write();
    hBkg_singleT_tW[iReg]->Write();
    hBkg_singleTbar_tW[iReg]->Write();
  
    out_Bkg[iReg]->Close();
  }
}
