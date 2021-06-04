void add_PUWeight_VtxCut(string MCname, string dataname, string MCNVtx, string dataNVtx,string PUweight){
  int nbins=65;
  float min=0;
  float max=65;
  TChain *tMC = new TChain("ntuplizer/tree");
  tMC->AddFile(MCname.c_str());
  tMC->AddFriend("tree",MCNVtx.c_str());
  TChain *tdata = new TChain("ntuplizer/tree");
  tdata->AddFile(dataname.c_str());
  tdata->AddFriend("tree",dataNVtx.c_str());
  string selection = "PV_N_good>0";
  vector<TH1F> hist_ratio;
  TH1F hist_ratio_1("hist_ratio_1","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_1);
  TH1F hist_ratio_2("hist_ratio_2","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_2);
  TH1F hist_ratio_3("hist_ratio_3","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_3);
  TH1F hist_ratio_4("hist_ratio_4","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_4);
  TH1F hist_ratio_5("hist_ratio_5","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_5);
  TH1F hist_ratio_6("hist_ratio_6","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_6);
  TH1F hist_ratio_7("hist_ratio_7","",nbins,min,max);
  hist_ratio.push_back(hist_ratio_7);

  vector<string> NPV_b {"PV_N_good","N_PV_inclu_cut","N_PV_verylow_cut","N_PV_low_cut","N_PV_medium_cut","N_PV_high_cut","N_PV_veryhigh_cut"};

  Float_t PUWeight[7];
  Float_t PUWeightE[7];
  TFile *f = new TFile(PUweight.c_str(),"recreate");
  TTree *T = new TTree("tree","PUweight");
  for(int i=0; i<7; i++){
    T->Branch(("PUWeight_"+NPV_b[i]).c_str(),&PUWeight[i],("PUWeight_"+NPV_b[i]+"/F").c_str());
    T->Branch(("PUWeightE_"+NPV_b[i]).c_str(),&PUWeightE[i],("PUWeightE_"+NPV_b[i]+"/F").c_str());
  }
  for(int ii=0; ii<7; ii++){
    TH1F *hist1 = new TH1F("hist1","Distribution of true pileup",nbins,min,max);
    TH1F *hist2 = new TH1F("hist2","Distribution of true pileup",nbins,min,max);
    TH1F *hist2_up = new TH1F("hist2_up","Distribution of true pileup",nbins,min,max);
    TH1F *hist2_down = new TH1F("hist2_down","Distribution of true pileup",nbins,min,max);
    TH1F *hist = new TH1F("hist","Distribution of true pileup",nbins,min,max);
    tMC->Draw("nPuVtxTrue>>hist",("("+selection+")*"+NPV_b[ii]).c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist1->Add(hist);
    tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selection+")*"+NPV_b[ii]).c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist2->Add(hist);

    tdata->Draw("Instan_Lumi_per_bunch_mean*69200*1.04>>hist",("("+selection+")*"+NPV_b[ii]).c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist2_up->Add(hist);
    tdata->Draw("Instan_Lumi_per_bunch_mean*69200*0.96>>hist",("("+selection+")*"+NPV_b[ii]).c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist2_down->Add(hist);

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
      hist_ratio[ii].SetBinContent(i,diff);
      hist_ratio[ii].SetBinError(i,diff_e);
    }
    delete hist1;
    delete hist2;
    delete hist2_up;
    delete hist2_down;
    delete hist;
  }

//  TFile * MCfile = new TFile(MCname.c_str(),"READ");
//  TTree * MCtree  = new TTree((TTree*)MCfile->Get("ntuplizer/tree"));

  vector<float> * nPuVtxTrue = new vector<float>[1];
//  nPuVtxTrue[0]= new float[5];
  nPuVtxTrue[0].push_back(0);
  tMC->SetBranchAddress("nPuVtxTrue",&nPuVtxTrue); 
  for (Int_t i=0; i<tMC->GetEntries(); i++){
    tMC->GetEntry(i);
    for(int ii=0; ii<7; ii++){
      PUWeight[ii] = hist_ratio[ii].GetBinContent(hist_ratio[ii].FindBin(nPuVtxTrue[0][0]));
      PUWeightE[ii] = hist_ratio[ii].GetBinError(hist_ratio[ii].FindBin(nPuVtxTrue[0][0]));
    }
    T->Fill();
  }
  T->Write();
  delete f;
}
