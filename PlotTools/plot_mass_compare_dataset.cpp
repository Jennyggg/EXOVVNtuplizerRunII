#include "helper.cpp"
void plot_mass_compare_dataset(string filenames_1,string filenames_2, string plotname, string legendname_1,string legendname_2,string massregion = "inclu", string insta_lumi = "inclu", bool isLog = false){
  int nbins=50;
  float min=0;
  float max=1000;
  std::vector<string> filename_1, filename_2;
//  splitstring(filenames_1, filename_1);
//  splitstring(filenames_2, filename_2);
  rootfilelist(filenames_1,filename_1);
  rootfilelist(filenames_2,filename_2);
  TChain *t = new TChain("ntuplizer/tree");
//  for(int i=0; i<filename_1.size(); i++){
//    t_1->AddFile(filename_1[i].c_str());
//  }
//  for(int i=0; i<filename_2.size(); i++){
//    t_2->AddFile(filename_2[i].c_str());
//  }
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
//  THStack *hs_1 = new THStack("hs_1","Distribution of mass/#tracks");
//  THStack *hs_2 = new THStack("hs_2","Distribution of mass/#tracks");
  THStack *hs_1 = new THStack("hs_1","Distribution of track mass");
  THStack *hs_2 = new THStack("hs_2","Distribution of track mass");
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  if(isLog){
    pad1->SetLogy();
//    hs_1->SetMinimum(0.1);
//    hs_2->SetMinimum(0.1);
  }
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  string selection = "!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0";
  if(massregion == "verylow") {selection += "&&Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40";max=40;min=20;nbins=20;}
  else if (massregion == "low") {selection += "&&Instanton_Trk_TrkCut_mass>40&&Instanton_Trk_TrkCut_mass<=80";max=80;min=40;nbins=40;}
  else if (massregion == "medium") {selection += "&&Instanton_Trk_TrkCut_mass>200&&Instanton_Trk_TrkCut_mass<=300";max=300;min=200;nbins=50;}
  else if (massregion == "high") {selection += "&&Instanton_Trk_TrkCut_mass>300&&Instanton_Trk_TrkCut_mass<=500";max=500;min=300;nbins=50;}
  else if (massregion != "inclu") {cout<<"Invalid mass region"<<endl; return;}

  if (insta_lumi != "inclu"){
    if(!(insta_lumi.at(2) == '-' &&isdigit(insta_lumi.at(0))&&isdigit(insta_lumi.at(1))&&isdigit(insta_lumi.at(3))&&isdigit(insta_lumi.at(4))))
      {cout<<"Invalid lumi range"<<endl; return;}
    selection += "&&Instan_Lumi_per_bunch_mean>0.000"+insta_lumi.substr(0,2)+"&&Instan_Lumi_per_bunch_mean<=0.000"+insta_lumi.substr(3,2);

  }
  TH1F *hist = new TH1F("hist","Distribution of track mass",nbins,min,max);
  TH1F *hist1 = new TH1F("hist1","Distribution of track mass",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of track mass",nbins,min,max);

  for(int i=0; i<filename_1.size(); i++){
    t->AddFile(filename_1[i].c_str());
    if(t->GetEntries(selection.c_str())>0){
      cout<<"processing file "<<filename_1[i]<<endl;
      t->Draw("Instanton_Trk_TrkCut_mass>>hist",selection.c_str(),"goff");
      hist=(TH1F*)gDirectory->Get("hist");
      hist1->Add(hist);
    }
    t->Reset();
  }
  for(int i=0; i<filename_2.size(); i++){
    t->AddFile(filename_2[i].c_str());
    if(t->GetEntries(selection.c_str())>0){
      cout<<"processing file "<<filename_2[i]<<endl;
      t->Draw("Instanton_Trk_TrkCut_mass>>hist",selection.c_str(),"goff");
      hist=(TH1F*)gDirectory->Get("hist");
      hist2->Add(hist);
    }
    t->Reset();
  }

//  t_1->Draw("(Instanton_Trk_TrkCut_mass)>>hist1",selection.c_str());
//  t_1->Draw("(Instanton_Trk_TrkCut_mass)>>hist1","Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&EVENT_run==316569&&Instan_Lumi_per_bunch_mean>0.00044&&Instan_Lumi_per_bunch_mean<=0.00045");
//  t_1->Draw("(Instanton_Trk_TrkCut_mass)>>hist1","Instanton_Trk_TrkCut_mass>0&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2");
//  t_1->Draw("(Instanton_Trk_TrkCut_mass)>>hist1","Instanton_Trk_TrkCut_mass>0&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&EVENT_run==316569&&Instan_Lumi_per_bunch_mean>0.00031&&Instan_Lumi_per_bunch_mean<=0.00032"); 
//  t_1->Draw("(Instanton_Trk_mass/Instanton_N_Trk)>>hist1","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); 
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","PV_N>20&&PV_N<25&&!Instanton_vtx_N_goodMuon&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

//  t_2->Draw("(Instanton_Trk_TrkCut_mass)>>hist2",selection.c_str());
//  t_2->Draw("(Instanton_Trk_TrkCut_mass)>>hist2","Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&EVENT_run==316569&&Instan_Lumi_per_bunch_mean>0.00044&&Instan_Lumi_per_bunch_mean<=0.00045");
//  t_2->Draw("(Instanton_Trk_TrkCut_mass)>>hist2","Instanton_Trk_TrkCut_mass>0&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2");
//  t_2->Draw("(Instanton_Trk_TrkCut_mass)>>hist2","Instanton_Trk_TrkCut_mass>0&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&EVENT_run==316569&&Instan_Lumi_per_bunch_mean>0.00031&&Instan_Lumi_per_bunch_mean<=0.00032");
//  t_2->Draw("(Instanton_Trk_mass/Instanton_N_Trk)>>hist2","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); 
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","PV_N>20&&PV_N<25&&!Instanton_vtx_N_goodMuon&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_2->Draw("Instanton_N_Trk_goodDisplaced>>hist2","!Instanton_vtx_isBPHtrigger&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total exclude BPH trigger collisions
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

//  hist1=(TH1F*)gDirectory->Get("hist1");
  hist1->Scale(1./hist1->GetSumOfWeights());
  hist1->SetFillColor(kOrange);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("Track Mass (GeV)");
  hs_1->GetYaxis()->SetTitle("Normalized Collisions");

//  hist2=(TH1F*)gDirectory->Get("hist2");
  hist2->Scale(1./hist2->GetSumOfWeights());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("Track Mass (GeV)");
  hs_2->GetYaxis()->SetTitle("Normalized Collisions");
  hs_1->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));
  hs_2->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist2,"UUNORM");

  TLegend *legend = new TLegend(0.65,0.65,0.90,0.90);
  legend->AddEntry(hist1,legendname_1.c_str(),"f");
  legend->AddEntry(hist2,legendname_2.c_str(),"lep");
  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();
  TLegend *legend1 = new TLegend(0.50,0.50,0.89,0.65);
  legend1->AddEntry((TObject*)0, "Instantaneous Luminosity per bunch:", "");
  if (insta_lumi == "inclu") legend1->AddEntry((TObject*)0, "   Inclusive Instantaneous Luminosity","");
  else legend1->AddEntry((TObject*)0, ("           0.000"+insta_lumi.substr(0,2)+"-0.000"+insta_lumi.substr(3,2)+" #mub^{-1}").c_str(),"");
//  legend1->AddEntry((TObject*)0, "           0.00044-0.00045 #mub^{-1}","");
//  legend1->AddEntry((TObject*)0, "   Inclusive Instantaneous Luminosity", "");
  if(massregion == "verylow")  legend1->AddEntry((TObject*)0, "           Very low mass region","");
  else if(massregion == "low")  legend1->AddEntry((TObject*)0, "           Low mass region","");
  else if(massregion == "medium")  legend1->AddEntry((TObject*)0, "           Medium mass region","");
  else if(massregion == "high")  legend1->AddEntry((TObject*)0, "           High mass region","");
  else legend1->AddEntry((TObject*)0, "           Inclusive mass region","");
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->Draw();

  pad2->cd();
/*  TH1F *hist_ratio = new TH1F("hist_ratio","",nbins,min,max);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->SetMarkerColor(kBlack);
  hist_ratio->GetXaxis()->SetTitle("Track Mass (GeV)");
  hist_ratio->GetYaxis()->SetTitle((legendname_2+"/"+legendname_1).c_str());
  hist_ratio->GetXaxis()->SetLabelSize(.06);
  hist_ratio->GetYaxis()->SetLabelSize(.06);
  hist_ratio->GetXaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleOffset(.55);
   for (Int_t i=0;i<nbins;i++) {
     Double_t diff;
     Double_t diff_e;
     if (hist1->GetBinContent(i)){
       diff = hist2->GetBinContent(i)/hist1->GetBinContent(i);
       diff_e = diff*sqrt(pow(hist1->GetBinError(i)/hist1->GetBinContent(i),2)+pow(hist2->GetBinError(i)/hist2->GetBinContent(i),2));
     }
     else{
       diff = 0; diff_e = 0;
     }
     hist_ratio->SetBinContent(i,diff);
     hist_ratio->SetBinError(i,diff_e);
   }
  hist_ratio->Draw("ep");
  hist_ratio->SetMinimum(0);
  hist_ratio->SetMaximum(std::max(2.,hist_ratio->GetMaximum()+hist_ratio->GetBinError(hist_ratio->GetMaximumBin())+0.1));
  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
*/
  TH1F *hist_diff = new TH1F("hist_diff","",nbins,min,max);
  hist_diff->SetMarkerStyle(21);
  hist_diff->SetMarkerColor(kBlack);
  hist_diff->GetXaxis()->SetTitle("Track Mass (GeV)");
  hist_diff->GetYaxis()->SetTitle("difference/#sigma(difference)");
  hist_diff->GetXaxis()->SetLabelSize(.06);
  hist_diff->GetYaxis()->SetLabelSize(.06);
  hist_diff->GetXaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleOffset(.55);
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
//    hist_diff->SetBinContent(i,diff/diff_e);
    hist_diff->SetBinError(i,0);
  }
  hist_diff->Draw("p");
//  hist_diff->SetMinimum(-5);
//  hist_ratio->SetMaximum(std::max(2.,hist_ratio->GetMaximum()+hist_ratio->GetBinError(hist_ratio->GetMaximumBin())+0.1));
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
  pad2->Update();
  pad1->Update();
  c->Update();
  plotname += "_"+massregion+"_mass_lumi"+insta_lumi;
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
