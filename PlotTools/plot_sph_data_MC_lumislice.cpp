#include "helper.cpp"
void plot_sph_data_MC_lumislice(string MCname, string dataname, string PUWeight_tree, string plotname, string MClegend, string datalegend, string MCPU="", string massregion = "inclu",string nPU_cut = "inclu", bool PUreweight = true, bool isLog = false, bool significanceplot = true){
  int nbins=20;
  float min=0.;
  float max=1.;
  TChain *tMC = new TChain("ntuplizer/tree");
  if(MCPU!="") tMC->AddFriend("tree",MCPU.c_str());
  TChain *tdata = new TChain("ntuplizer/tree");
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  THStack *hs_1 = new THStack("hs_1","Distribution of Spherocity");
  THStack *hs_2 = new THStack("hs_2","Distribution of Spherocity");
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
  string weight;
  string weighte;
  string selection = "!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0";
  if(massregion == "verylow") {selection += "&&Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&Instanton_Trk_TrkCut_spherocity<0.3";weight="PUWeight_N_PV_verylow_cut";weighte="PUWeightE_N_PV_verylow_cut";}
  else if (massregion == "low") {selection += "&&Instanton_Trk_TrkCut_mass>40&&Instanton_Trk_TrkCut_mass<=80&&Instanton_Trk_TrkCut_spherocity<0.3";weight="PUWeight_N_PV_low_cut";weighte="PUWeightE_N_PV_low_cut";}
  else if (massregion == "medium") {selection += "&&Instanton_Trk_TrkCut_mass>80&&Instanton_Trk_TrkCut_mass<=200";weight="PUWeight_N_PV_medium_cut";weighte="PUWeightE_N_PV_medium_cut";}
  else if (massregion == "high") {selection += "&&Instanton_Trk_TrkCut_mass>200&&Instanton_Trk_TrkCut_mass<=300";weight="PUWeight_N_PV_high_cut";weighte="PUWeightE_N_PV_high_cut";}
  else if (massregion == "veryhigh") {selection += "&&Instanton_Trk_TrkCut_mass>300&&Instanton_Trk_TrkCut_mass<=500";weight="PUWeight_N_PV_veryhigh_cut";weighte="PUWeightE_N_PV_veryhigh_cut";}
  else if (massregion != "inclu") {cout<<"Invalid mass region"<<endl; return;}
  weight="PUWeight_N_PV_inclu_cut";weighte="PUWeightE_N_PV_inclu_cut";
  string selectionMC,selectiondata;
  selectionMC=selection;
  selectiondata=selection;

  if (nPU_cut != "inclu"){
    if(!(nPU_cut.at(2) == '-' &&(isdigit(nPU_cut.at(0))||nPU_cut.at(0)== ' ')&&(isdigit(nPU_cut.at(1))||nPU_cut.at(0)==' ')&&isdigit(nPU_cut.at(3))&&isdigit(nPU_cut.at(4))))
      {cout<<"Invalid lumi range"<<endl; return;}
    selectionMC += "&&nPuVtxTrue_>"+nPU_cut.substr(0,2)+"&&nPuVtxTrue_<="+nPU_cut.substr(3,2);
    selectiondata += "&&Instan_Lumi_per_bunch_mean*69200>"+nPU_cut.substr(0,2)+"&&Instan_Lumi_per_bunch_mean*69200<="+nPU_cut.substr(3,2);
  }
  cout<<"selectionMC "<<selectionMC<<endl;
  cout<<"selectiondata "<<selectiondata<<endl;
  TH1F *hist = new TH1F("hist","Distribution of Spherocity",nbins,min,max);
  TH1F *hist1 = new TH1F("hist1","Distribution of Spherocity",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of Spherocity",nbins,min,max);
  hist->Sumw2(kTRUE);
  tMC->AddFile(MCname.c_str());
  tMC->AddFriend("tree",PUWeight_tree.c_str());
  if(PUreweight){
    if(tMC->GetBranchStatus("Instanton_vtx_genvertex_PVDistance"))
      tMC->Draw("Instanton_Trk_TrkCut_spherocity>>hist",("("+selectionMC+"&&Instanton_vtx_genvertex_PVDistance>2)*"+weight).c_str(),"goff");
    else
      tMC->Draw("Instanton_Trk_TrkCut_spherocity>>hist",("("+selectionMC+")*"+weight).c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist1->Add(hist);
    if(tMC->GetBranchStatus("Instanton_vtx_genvertex_PVDistance"))
      HistError(tMC,PUWeight_tree,"Instanton_Trk_TrkCut_spherocity",weighte,selectionMC+"&&Instanton_vtx_genvertex_PVDistance>2",hist1,nbins,min,max);
    else
      HistError(tMC,PUWeight_tree,"Instanton_Trk_TrkCut_spherocity",weighte,selectionMC,hist1,nbins,min,max);
  }
  else{
    if(tMC->GetBranchStatus("Instanton_vtx_genvertex_PVDistance"))
      tMC->Draw("Instanton_Trk_TrkCut_spherocity>>hist",("("+selectionMC+"&&Instanton_vtx_genvertex_PVDistance>2)").c_str(),"goff");
    else
      tMC->Draw("Instanton_Trk_TrkCut_spherocity>>hist",("("+selectionMC+")").c_str(),"goff");
    hist=(TH1F*)gDirectory->Get("hist");
    hist1->Add(hist);
  }

  tdata->AddFile(dataname.c_str());
  tdata->Draw("Instanton_Trk_TrkCut_spherocity>>hist",selectiondata.c_str(),"goff");
  hist=(TH1F*)gDirectory->Get("hist");
  hist2->Add(hist);

  for(int i=0; i<=nbins; i++){
    cout<<"bin "<<i<<", hist1 unnormalize "<<hist1->GetBinContent(i)<<", error "<<hist1->GetBinError(i)<<endl;
  }
  cout<<"hist1 SumWeight "<<hist1->GetSumOfWeights()<<endl;
  cout<<"hist1 Entries "<<hist1->GetEntries()<<endl;

  hist1->Scale(1./hist1->GetSumOfWeights());
  hist1->SetFillColor(kOrange);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);

  cout<<"nPuVtxTrue_ "<<(int)tMC->GetBranchStatus("nPuVtxTrue_")<<", max = "<<tMC->GetMaximum("nPuVtxTrue_")<<", min = "<<tMC->GetMinimum("nPuVtxTrue_")<<endl;
  cout<<"Instan_Lumi_per_bunch_mean "<<(int)tdata->GetBranchStatus("Instan_Lumi_per_bunch_mean")<<", max = "<<tdata->GetMaximum("Instan_Lumi_per_bunch_mean")*69200<<", min = "<<tdata->GetMinimum("Instan_Lumi_per_bunch_mean")*69200<<endl;

  for(int i=0; i<=nbins; i++){
    cout<<"bin "<<i<<", hist1 "<<hist1->GetBinContent(i)<<", error "<<hist1->GetBinError(i)<<endl;
  }
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("Spherocity");
  hs_1->GetYaxis()->SetTitle("Normalized Collisions");

  for(int i=0; i<=nbins; i++){
    cout<<"bin "<<i<<", hist2 unnormalize "<<hist2->GetBinContent(i)<<", error "<<hist2->GetBinError(i)<<endl;
  }
  cout<<"hist2 SumWeight "<<hist2->GetSumOfWeights()<<endl;
  cout<<"hist2 Entries "<<hist2->GetEntries()<<endl;
  hist2->Scale(1./hist2->GetSumOfWeights());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  for(int i=0; i<=nbins; i++){
    cout<<"bin "<<i<<", hist2 "<<hist2->GetBinContent(i)<<", error "<<hist2->GetBinError(i)<<endl;
  }
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("Spherocity");
  hs_2->GetYaxis()->SetTitle("Normalized Collisions");
  hs_1->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));
  hs_2->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist2,"UUNORM");

  TLegend *legend = new TLegend(0.65,0.65,0.90,0.90);
  legend->AddEntry(hist1,MClegend.c_str(),"f");
  legend->AddEntry(hist2,datalegend.c_str(),"lep");
  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();
  TLegend *legend1 = new TLegend(0.50,0.50,0.89,0.65);
//  legend1->AddEntry((TObject*)0, "Instantaneous Luminosity per bunch:", "");
  if (nPU_cut == "inclu") legend1->AddEntry((TObject*)0, "   Inclusive Pileup","");
  else legend1->AddEntry((TObject*)0, ("           N truePU   "+nPU_cut.substr(0,2)+"-"+nPU_cut.substr(3,2)).c_str(),"");
  if(massregion == "verylow")  legend1->AddEntry((TObject*)0, "           Very low mass region","");
  else if(massregion == "low")  legend1->AddEntry((TObject*)0, "           Low mass region","");
  else if(massregion == "medium")  legend1->AddEntry((TObject*)0, "           Medium mass region","");
  else if(massregion == "high")  legend1->AddEntry((TObject*)0, "           High mass region","");
  else if(massregion == "veryhigh")  legend1->AddEntry((TObject*)0, "         Very High mass region","");
  else legend1->AddEntry((TObject*)0, "           Inclusive mass region","");
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->Draw();

  pad2->cd();
  TH1F *hist_diff = new TH1F("hist_diff","",nbins,min,max);
  hist_diff->SetMarkerStyle(21);
  hist_diff->SetMarkerColor(kBlack);
  hist_diff->GetXaxis()->SetTitle("Spherocity");
  if(significanceplot)
    hist_diff->GetYaxis()->SetTitle("difference/#sigma(difference)");
  else
    hist_diff->GetYaxis()->SetTitle("ratio");
  hist_diff->GetXaxis()->SetLabelSize(.06);
  hist_diff->GetYaxis()->SetLabelSize(.06);
  hist_diff->GetXaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleOffset(.55);
  if(significanceplot){
    for (Int_t i=0;i<nbins;i++) {
      Double_t diff;
      Double_t diff_e;
      if (hist1->GetBinContent(i)||hist2->GetBinContent(i)){
        diff = hist2->GetBinContent(i)-hist1->GetBinContent(i);
        diff_e = sqrt(hist1->GetBinError(i)*hist1->GetBinError(i)+hist2->GetBinError(i)*hist2->GetBinError(i));
        hist_diff->SetBinContent(i,diff/diff_e);
      }
      else{
        diff = 0; diff_e = 0;
        hist_diff->SetBinContent(i,0);
      }
      hist_diff->SetBinError(i,0);
    }
  }
  else{
    for (Int_t i=0;i<nbins;i++) {
      Double_t diff;
      Double_t diff_e;
      if (hist1->GetBinContent(i)&&hist2->GetBinContent(i)){
        diff = hist2->GetBinContent(i)/hist1->GetBinContent(i);
        diff_e = diff*sqrt(pow(hist1->GetBinError(i)/hist1->GetBinContent(i),2)+pow(hist2->GetBinError(i)/hist2->GetBinContent(i),2));
      }
      else{
        diff = 0; diff_e = 0;
      }
      hist_diff->SetBinContent(i,diff);
      hist_diff->SetBinError(i,diff_e);
    }
  }
  hist_diff->Draw("p");
  if(significanceplot){
    TLine *line = new TLine(min,0,max,0);
    line->SetLineColor(kRed);
    line->Draw("SAME");
    TLine *line_up = new TLine(min,1,max,1);
    line_up ->SetLineStyle(9);
    line_up->SetLineColor(kBlue);
    TLine *line_down = new TLine(min,-1,max,-1);
    line_down ->SetLineStyle(9);
    line_down->SetLineColor(kBlue);
    line_up->Draw("SAME");
    line_down->Draw("SAME");
  }
  else{
    TLine *line = new TLine(min,1,max,1);
    line->SetLineColor(kRed);
    line->Draw("SAME");
  }
  pad2->Update();
  pad1->Update();
  c->Update();
  plotname += "_"+massregion+"_PU"+nPU_cut+"_sph";
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
  delete c;

}
