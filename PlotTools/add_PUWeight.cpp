void add_PUWeight(string MCname, string dataname,string PUweight){
  int nbins=65;
  float min=0;
  float max=65;
  TH1F *hist1 = new TH1F("hist1","Distribution of true pileup",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of true pileup",nbins,min,max);
  TH1F *hist2_up = new TH1F("hist2_up","Distribution of true pileup",nbins,min,max);
  TH1F *hist2_down = new TH1F("hist2_down","Distribution of true pileup",nbins,min,max);
  TH1F *hist = new TH1F("hist","Distribution of true pileup",nbins,min,max);
  TChain *tMC = new TChain("ntuplizer/tree");
  TChain *tdata = new TChain("ntuplizer/tree");
  string selection = "PV_N_good>0";

  tMC->AddFile(MCname.c_str());
  tMC->Draw("nPuVtxTrue>>hist",selection.c_str(),"goff");
  hist=(TH1F*)gDirectory->Get("hist");
  hist1->Add(hist);

  tdata->AddFile(dataname.c_str());
  tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",selection.c_str(),"goff");
  hist=(TH1F*)gDirectory->Get("hist");
  hist2->Add(hist);

  tdata->Draw("Instan_Lumi_per_bunch_mean*69200*1.04>>hist",selection.c_str(),"goff");
  hist=(TH1F*)gDirectory->Get("hist");
  hist2_up->Add(hist);
  tdata->Draw("Instan_Lumi_per_bunch_mean*69200*0.96>>hist",selection.c_str(),"goff");
  hist=(TH1F*)gDirectory->Get("hist");
  hist2_down->Add(hist);

  TH1F *hist_ratio = new TH1F("hist_ratio","",nbins,min,max);
  for (Int_t i=0;i<=nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    Double_t diff_up;
    Double_t diff_down;
    Double_t lumi_e;
    if (hist1->GetBinContent(i)&&hist2->GetBinContent(i)){
      diff = hist2->GetBinContent(i)/hist1->GetBinContent(i);
      diff_e = diff*sqrt(pow(hist1->GetBinError(i)/hist1->GetBinContent(i),2)+pow(hist2->GetBinError(i)/hist2->GetBinContent(i),2));
      diff_up = hist2_up->GetBinContent(i)/hist1->GetBinContent(i);
      diff_down = hist2_down->GetBinContent(i)/hist1->GetBinContent(i);
      lumi_e = std::max(fabs(diff_up-diff),fabs(diff_down-diff));
      diff_e = sqrt(diff_e*diff_e+lumi_e*lumi_e);
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist_ratio->SetBinContent(i,diff);
    hist_ratio->SetBinError(i,diff_e);
  }

  Float_t PUWeight;
  Float_t PUWeightE;
  TFile *f = new TFile(PUweight.c_str(),"recreate");
  TTree *T = new TTree("tree","PUweight");
  T->Branch("PUWeight",&PUWeight,"PUWeight/F");
  T->Branch("PUWeightE",&PUWeightE,"PUWeightE/F");

//  TFile * MCfile = new TFile(MCname.c_str(),"READ");
//  TTree * MCtree  = new TTree((TTree*)MCfile->Get("ntuplizer/tree"));

  vector<float> * nPuVtxTrue = new vector<float>[1];
//  nPuVtxTrue[0]= new float[5];
  nPuVtxTrue[0].push_back(0);
  tMC->SetBranchAddress("nPuVtxTrue",&nPuVtxTrue); 
  for (Int_t i=0; i<tMC->GetEntries(); i++){
    tMC->GetEntry(i);
    PUWeight = hist_ratio->GetBinContent(hist_ratio->FindBin(nPuVtxTrue[0][0]));
    PUWeightE = hist_ratio->GetBinError(hist_ratio->FindBin(nPuVtxTrue[0][0]));
    T->Fill();
  }
  T->Write();
  delete f;
}
