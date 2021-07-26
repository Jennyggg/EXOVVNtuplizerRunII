#include "helper.cpp"


string OBS_default="thrust";
string MASSREGION_default="inclu";
//string SIGNALNAME_default="../Ntuplizer/flatTuple_h7_4f_noPU.root ../Ntuplizer/flatTuple_h7_5f_noPU.root ../Ntuplizer/flatTuple_h7_4f.root ../Ntuplizer/flatTuple_h7_5f.root";
//string SIGNALNAME_default = "../Ntuplizer/flatTuple_sherpa_4f.root ../Ntuplizer/flatTuple_sherpa_5f.root";
string SIGNALNAME_default = "../Ntuplizer/flatTuple_h7_4f.root ../Ntuplizer/flatTuple_sherpa_4f.root";
//string SIGNALNAME_default = "../Ntuplizer/flatTuple_h7_5f.root ../Ntuplizer/flatTuple_sherpa_5f.root";
string BKGNAME_default="/pnfs/psi.ch/cms/trivcat/store/user/jinw/MinBiasMC/CRAB/MinBias_NoFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8/MinBias_NoFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/210622_211735/0000/flatTuple.root";
//string PLOTNAME_default="plots_sig_bkg/Sig_sherpa_Bkg";
//string PLOTNAME_default="plots_sig_bkg/Sig_h7PUcompare_Bkg";
string PLOTNAME_default="plots_sig_bkg/Sig_4f_h7_sherpa_Bkg";
//string PLOTNAME_default="plots_sig_bkg/Sig_5f_h7_sherpa_Bkg";
//string SIGNALLEGEND_default="Instanton-sherpa-4f-withPU Instanton-sherpa-5f-withPU";
//string SIGNALLEGEND_default="Instanton-h7-4f-noPU Instanton-h7-5f-noPU Instanton-h7-4f-withPU Instanton-h7-5f-withPU";
string SIGNALLEGEND_default="Instanton-h7-4f-withPU Instanton-sherpa-4f-withPU";
string BKGLEGEND_default="MinBias";

map<string,pair<double,double>> XS{
{"Herwig-rand3-fix4f",{5.81406e+09, 87884.9e+03}},
{"Herwig-rand3-fix5f",{5.76513e+09, 86810.3e+03}},
{"Sherpa-fix4f-scale1",{4.97368e+10,4.2934e+08}},
{"Sherpa-fix4f-scale2",{4.18972e+09,3.61958e+07}},
{"Sherpa-fix4f-scale0p5",{2.35703e+08,1.15072e+06}},
{"Sherpa-fix5f-scale1",{4.90225e+10,4.2277e+08}},
{"Sherpa-fix5f-scale2",{7.47231e+10,6.36042e+08}},
{"Sherpa-fix5f-scale0p5",{2.30001e+08,1.47157e+06}},
{"MinBias",{69.2e+09,3.1832e+09}}
};
vector<double> SIGNALXS_h7PUcompare = {XS["Herwig-rand3-fix4f"].first,XS["Herwig-rand3-fix5f"].first,XS["Herwig-rand3-fix4f"].first,XS["Herwig-rand3-fix5f"].first};
vector<double> SIGNALXS_E_h7PUcompare = {XS["Herwig-rand3-fix4f"].second,XS["Herwig-rand3-fix5f"].second,XS["Herwig-rand3-fix4f"].second,XS["Herwig-rand3-fix5f"].second};
vector<double> SIGNALXS_4f_h7_sherpa = {XS["Herwig-rand3-fix4f"].first, XS["Sherpa-fix4f-scale1"].first};
vector<double> SIGNALXS_E_4f_h7_sherpa = {XS["Herwig-rand3-fix4f"].second, XS["Sherpa-fix4f-scale1"].second};
vector<double> SIGNALXS_5f_h7_sherpa = {XS["Herwig-rand3-fix5f"].first, XS["Sherpa-fix5f-scale1"].first};
vector<double> SIGNALXS_E_5f_h7_sherpa = {XS["Herwig-rand3-fix5f"].second, XS["Sherpa-fix5f-scale1"].second};
vector<double> SIGNALXS_sherpa = {XS["Sherpa-fix4f-scale1"].first, XS["Sherpa-fix5f-scale1"].first};
vector<double> SIGNALXS_E_sherpa = {XS["Sherpa-fix4f-scale1"].second, XS["Sherpa-fix5f-scale1"].second};

vector<double> SIGNALXS_default = SIGNALXS_4f_h7_sherpa;
vector<double> SIGNALXS_E_default = SIGNALXS_E_4f_h7_sherpa;
vector<double> BKGXS_default = {XS["MinBias"].first};
vector<double> BKGXS_E_default = {XS["MinBias"].second};
bool ISLOG_default = true;



TH1F* MakeHist_Sig(TChain * tree, string selection,string histname, string Obs, double max, double min, int nbins, double XS, double e_XS, Color_t color, TTree * friendtree=0){
  if(friendtree!=0) tree->AddFriend(friendtree);
  TH1F *hist = new TH1F(histname.c_str(),histname.c_str(),nbins,min,max);
  hist->Sumw2(kTRUE);
  tree->Draw((BRANCHNAMEMC[Obs]+">>"+histname).c_str(),("("+selection+")&&genWeight").c_str(),"goff");
//  tree->Draw((BRANCHNAMEMC[Obs]+">>"+histname).c_str(),("("+selection+")").c_str(),"goff");
  hist=(TH1F*)gDirectory->Get(histname.c_str());
  if(hist->GetSumOfWeights()==0){
    cout<<"The sum of weight is zero !"<<endl;
    return 0;
  }
  double ScaleF = XS/hist->GetSumOfWeights();
//  cout<<"sig"<<endl;
  for(int i=0; i<=nbins; i++){
    if(hist->GetBinContent(i)>0 && XS>0){
//      cout<<"before "<<hist->GetBinContent(i)<<" "<<hist->GetBinError(i)<<" "<<XS<<" "<<e_XS<<endl;
      hist->SetBinError(i,ScaleF*hist->GetBinContent(i)*sqrt(e_XS*e_XS/XS/XS+hist->GetBinError(i)*hist->GetBinError(i)/hist->GetBinContent(i)/hist->GetBinContent(i)));
      hist->SetBinContent(i,hist->GetBinContent(i)*ScaleF);
//      cout<<"after "<<hist->GetBinContent(i)<<" "<<hist->GetBinError(i)<<" "<<XS<<" "<<e_XS<<endl;
    }
  }
  hist->SetLineColor(color);
  hist->SetMarkerStyle(21);
  hist->SetMarkerColor(color);
  hist->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist->GetYaxis()->SetTitle("Normalized Collisions");
  hist->Draw("elpSAME");
  return hist;
}


TH1F* MakeHist_Bkg(TChain * tree, string selection,string histname, string Obs, double max, double min, int nbins, double XS, double e_XS, Color_t color, TTree * friendtree=0){
  if(friendtree!=0) tree->AddFriend(friendtree);
  TH1F *hist = new TH1F(histname.c_str(),histname.c_str(),nbins,min,max);
  hist->Sumw2(kTRUE);
  tree->Draw((BRANCHNAMEMC[Obs]+">>"+histname).c_str(),("("+selection+")&&genWeight").c_str(),"goff");
  hist=(TH1F*)gDirectory->Get(histname.c_str());
  if(hist->GetSumOfWeights()==0){
    cout<<"The sum of weight is zero !"<<endl;
    return 0;
  }
  double ScaleF = XS/hist->GetSumOfWeights();
//  cout<<"bkg"<<endl;
  for(int i=0; i<=nbins; i++){
    if(hist->GetBinContent(i)>0 && XS>0){
//      cout<<"before "<<hist->GetBinContent(i)<<" "<<hist->GetBinError(i)<<" "<<XS<<" "<<e_XS<<endl;
      hist->SetBinError(i,ScaleF*hist->GetBinContent(i)*sqrt(e_XS*e_XS/XS/XS+hist->GetBinError(i)*hist->GetBinError(i)/hist->GetBinContent(i)/hist->GetBinContent(i)));
      hist->SetBinContent(i,hist->GetBinContent(i)*ScaleF);
//      cout<<"after "<<hist->GetBinContent(i)<<" "<<hist->GetBinError(i)<<" "<<XS<<" "<<e_XS<<endl;
    }
  }
  hist->SetFillColor(color);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist->GetYaxis()->SetTitle("Normalized Collisions");
//  hist->Draw("HIST");
//  hist->Draw("ESAME");
  return hist;
}

TH1F* MakeHist_Ratio(TH1F* SigHist, TH1F * BkgHist, string histname,string Obs, double max, double min, int nbins){
  TH1F *hist = new TH1F(histname.c_str(),"",nbins,min,max);
  hist->SetMarkerStyle(21);
  hist->SetMarkerColor(SigHist->GetMarkerColor());
  hist->SetLineColor(SigHist->GetLineColor());
  hist->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist->GetYaxis()->SetTitle("ratio");
  hist->GetXaxis()->SetLabelSize(.06);
  hist->GetYaxis()->SetLabelSize(.06);
  hist->GetXaxis()->SetTitleSize(.08);
  hist->GetYaxis()->SetTitleSize(.08);
  hist->GetYaxis()->SetTitleOffset(.55);
  for (Int_t i=0;i<nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    if (SigHist->GetBinContent(i)&&BkgHist->GetBinContent(i)){
      diff = SigHist->GetBinContent(i)/BkgHist->GetBinContent(i);
      diff_e = diff*sqrt(pow(BkgHist->GetBinError(i)/BkgHist->GetBinContent(i),2)+pow(SigHist->GetBinError(i)/SigHist->GetBinContent(i),2));
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist->SetBinContent(i,diff);
    hist->SetBinError(i,diff_e);
  }
  hist->SetMaximum(2.);
  hist->SetMinimum(0.);
  hist->Draw("elpSAME");
  return hist;
}

void plot_obs_sigs_bkgs(string Obs=OBS_default,string massregion=MASSREGION_default,string SignalName=SIGNALNAME_default, string BkgName=BKGNAME_default, string plotname=PLOTNAME_default,string SignalLegend=SIGNALLEGEND_default,string BkgLegend=BKGLEGEND_default,vector<double> SignalXS = SIGNALXS_default,vector<double> SignalXS_E=SIGNALXS_E_default,vector<double> BkgXS = BKGXS_default,vector<double> BkgXS_E=BKGXS_E_default,bool isLog = ISLOG_default){

  vector<string> SignalNames;
  vector<string> BkgNames;
  vector<string> SignalLegends;
  vector<string> BkgLegends;
  splitstring(SignalName,SignalNames);
  splitstring(BkgName,BkgNames);
  splitstring(SignalLegend, SignalLegends);
  splitstring(BkgLegend, BkgLegends);
  vector<Color_t> SigColors={kRed,kGreen,kViolet,kCyan,kMagenta,kYellow};
  vector<Color_t> BkgColors={kOrange,kSpring+2,kCyan-6, kGray};
  if(Obs!="mass"&&Obs!="pu"&&Obs!="sph"&&Obs!="thrust"&&Obs!="ndisplace"&&Obs!="ndisplaceQ4"&&Obs!="ntrk"&&Obs!="ntrk_pt05"&&Obs!="ntrk_pt1"&&Obs!="pt"&&Obs!="pt_pt05"&&Obs!="pt_pt1"&&Obs!="masspertrk"&&Obs!="nSVsig2"&&Obs!="nSVsig3"&&Obs!="nSVsig5"&&Obs!="ntrkSVsig2"&&Obs!="ntrkSVsig3"&&Obs!="ntrkSVsig5"&&Obs!="ntrkdisplaceV"&&Obs!="ntrkdisplaceV_svsig2"){
    cout<<"Invalid Observable"<<endl;
    return;
  }
  if(massregion!="incluCR"&&massregion!="inclu"&&massregion!="verylow"&&massregion!="verylowfull"&&massregion!="low"&&massregion!="lowfull"&&massregion!="medium"&&massregion!="high"&&massregion!="veryhigh"){
    cout<<"Invalid mass region"<<endl;
    return;
  }
  int nbins=NBINS[Obs];
  float min=MIN[Obs];
  float max=MAX[Obs];
  vector<TChain *> tSignals;
  vector<TTree *> tSignalFriends;
  for(int isig=0; isig < SignalNames.size(); isig++){
    tSignals.push_back(new TChain("ntuplizer/tree"));
    tSignals[isig]->AddFile(SignalNames[isig].c_str());
    tSignalFriends.push_back(new TTree("tree","isClosestPV"));
    vector<bool> isClosestPV_signal;
    tSignalFriends[isig]->Branch("isClosestPV",&isClosestPV_signal);
    vector<float> *PV_x_signal=0;
    vector<float> *PV_y_signal=0;
    vector<float> *PV_z_signal=0;
    float ins_VX,ins_VY,ins_VZ;
    tSignals[isig]->SetBranchAddress("PV_x",&PV_x_signal);
    tSignals[isig]->SetBranchAddress("PV_y",&PV_y_signal);
    tSignals[isig]->SetBranchAddress("PV_z",&PV_z_signal);
    tSignals[isig]->SetBranchAddress("Instanton_VX",&ins_VX);
    tSignals[isig]->SetBranchAddress("Instanton_VY",&ins_VY);
    tSignals[isig]->SetBranchAddress("Instanton_VZ",&ins_VZ);
//Find recovertex closest to the genvertex and label it.
    for(int i=0; i<tSignals[isig]->GetEntries(); i++){
      float min_distance2=999.;
      float distance2=0;
      int min_distance2_index= -1;
      tSignals[isig]->GetEntry(i);
      isClosestPV_signal.insert(isClosestPV_signal.end(),PV_x_signal->size(),false);
      for(size_t ii=0; ii<PV_x_signal->size(); ii++){
        distance2=(PV_x_signal->at(ii)-ins_VX)*(PV_x_signal->at(ii)-ins_VX)+(PV_y_signal->at(ii)-ins_VY)*(PV_y_signal->at(ii)-ins_VY)+(PV_z_signal->at(ii)-ins_VZ)*(PV_z_signal->at(ii)-ins_VZ);
        if(distance2<min_distance2){
          min_distance2=distance2;
          min_distance2_index=(int)ii;
        }
      }
      isClosestPV_signal[min_distance2_index]=true;
      tSignalFriends[isig]->Fill();
      isClosestPV_signal.clear();
    }
  }

  vector<TChain *> tBkgs;
  for(int ibkg=0; ibkg < BkgNames.size(); ibkg++){
    tBkgs.push_back(new TChain("ntuplizer/tree"));
    tBkgs[ibkg]->AddFile(BkgNames[ibkg].c_str());
  }
  string selection = SELE_MASS[massregion];
  string selectionSignal,selectionBkg;
// Select the vertex that satify the selection of mass and is closest to the genvertex for each event.
  selectionSignal="("+selection+")&&isClosestPV";
// For the background all the vertices including pileup are selected if they satisfy the selestions for the mass regions.
  selectionBkg=selection;
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  if(isLog){
    pad1->SetLogy();
  }
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
// Draw backgrouds in stacked historgrams (for now we only have MinBias background).
  THStack *hs_Bkg = new THStack("hs_Bkg",PLOTNAME[Obs].c_str());
  vector<TH1F *> hs_Bkgs;
  for(int ibkg=0; ibkg<tBkgs.size(); ibkg++){
    hs_Bkgs.push_back(MakeHist_Bkg(tBkgs[ibkg], selectionBkg,"Bkg"+to_string(ibkg), Obs, max, min, nbins, BkgXS[ibkg], BkgXS_E[ibkg], BkgColors[ibkg]));
    hs_Bkg->Add(hs_Bkgs[ibkg]);
  }
  hs_Bkg->Draw("hist");
  hs_Bkg->Draw("ESAME");
  hs_Bkg->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_Bkg->GetYaxis()->SetTitle("Collisions normalized to 1 pb^{-1}");
// Draw signals.
  vector<TH1F *> hs_Sigs;
  for(int isig=0; isig<tSignals.size(); isig++){
    hs_Sigs.push_back(MakeHist_Sig(tSignals[isig], selectionSignal,"Sig"+to_string(isig), Obs, max, min, nbins, SignalXS[isig], SignalXS_E[isig], SigColors[isig],tSignalFriends[isig]));
  }
//Change the min and max of the plots according to signal and background range.
  vector<double> sig_bkg_min;
  for(int isig=0; isig<tSignals.size(); isig++){
    sig_bkg_min.push_back(0.1*hs_Sigs[isig]->GetSumOfWeights()/hs_Sigs[isig]->GetEntries());
  }
  hs_Bkg->SetMinimum(*min_element(sig_bkg_min.begin(),sig_bkg_min.end()));

  vector<double> sig_bkg_max;
  for(int isig=0; isig<tSignals.size(); isig++){
    sig_bkg_max.push_back(hs_Sigs[isig]->GetBinContent(hs_Sigs[isig]->GetMaximumBin()));
  }
  sig_bkg_max.push_back(((TH1F*)hs_Bkg->GetStack()->Last())->GetBinContent(((TH1F*)hs_Bkg->GetStack()->Last())->GetMaximumBin()));
  hs_Bkg->SetMaximum(std::max((double)((TH1F*)hs_Bkg->GetStack()->Last())->GetMaximum(),*max_element(sig_bkg_max.begin(),sig_bkg_max.end())*1.5));



  TLegend *legend = new TLegend(0.65,0.65,0.90,0.90);
  for(int isig=0; isig<tSignals.size(); isig++){
    legend->AddEntry(hs_Sigs[isig],SignalLegends[isig].c_str(),"lep");
  }
  for(int ibkg=0; ibkg<tBkgs.size(); ibkg++){
    legend->AddEntry(hs_Bkgs[ibkg],BkgLegends[ibkg].c_str(),"f");
  }
  legend->SetFillStyle(0);
  legend->Draw();
  TLegend *legend1 = new TLegend(0.50,0.50,0.89,0.65);
  legend1->AddEntry((TObject*)0,LEGENDNAME[massregion].c_str(),"");
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->Draw();
  pad1->Update();

  pad2->cd();
// Plot the signal to background ratio
  vector<TH1F*> hs_Ratios;
  for(int isig=0; isig<tSignals.size(); isig++){
    hs_Ratios.push_back(MakeHist_Ratio(hs_Sigs[isig], (TH1F*)hs_Bkg->GetStack()->Last(), "s/b"+to_string(isig),Obs, max, min, nbins));
  }
  vector<double> sig_bkg_ratio_max;
  for(int isig=0; isig<tSignals.size(); isig++){
    sig_bkg_ratio_max.push_back(hs_Ratios[isig]->GetBinContent(hs_Ratios[isig]->GetMaximumBin()));
  }
  hs_Ratios[0]->SetMaximum(std::max(2.,*max_element(sig_bkg_ratio_max.begin(),sig_bkg_ratio_max.end())));

  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kBlue);
  line->Draw("SAME");
  pad2->Update();
  c->Update();
  plotname += "_"+massregion+"_PUinclu"+"_"+Obs;
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());


}
