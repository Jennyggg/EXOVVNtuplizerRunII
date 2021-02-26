#include "helper.cpp"
void plot_ndisplace_compare_trig(string filenames_1,string filenames_2, string plotname, string legendname_1,string legendname_2, bool isLog = false){
  int nbins=20;
  float min=0.;
  float max=20;
  std::vector<string> filename_1, filename_2;
  splitstring(filenames_1, filename_1);
  splitstring(filenames_2, filename_2);
  TChain *t_1 = new TChain("ntuplizer/tree");
  TChain *t_2 = new TChain("ntuplizer/tree");
  for(int i=0; i<filename_1.size(); i++){
    t_1->AddFile(filename_1[i].c_str());
  }
  for(int i=0; i<filename_2.size(); i++){
    t_2->AddFile(filename_2[i].c_str());
  }
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  THStack *hs_1 = new THStack("hs_1","Distribution of #displaced tracks");
  THStack *hs_2 = new THStack("hs_2","Distribution of #displaced tracks");
  THStack *hs_3 = new THStack("hs_3","Distribution of #displaced tracks");
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
  TH1F *hist1 = new TH1F("hist1","Distribution of #displaced tracks",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of #displaced tracks",nbins,min,max);
  TH1F *hist3 = new TH1F("hist3","Distribution of #displaced tracks",nbins,min,max);

  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","!Instanton_vtx_isBPHtrigger&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total exclude BPH trigger collisions
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

  t_2->Draw("Instanton_N_Trk_Displaced>>hist3","Instanton_vtx_isBPHtrigger&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total BPH trigger collisions

  hist1=(TH1F*)gDirectory->Get("hist1");
  hist1->Scale(1./hist1->GetEntries());
  hist1->SetFillColor(kOrange);
//  hist1->SetLineColor(kBlack);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("# Displaced Tracks");
  hs_1->GetYaxis()->SetTitle("Normalized Collisions");

  hist2=(TH1F*)gDirectory->Get("hist2");
  hist2->Scale(1./hist2->GetEntries());
  hist2->SetLineColor(kBlue);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("# Displaced Tracks");
  hs_2->GetYaxis()->SetTitle("Normalized Collisions");

  hist3=(TH1F*)gDirectory->Get("hist3");
  hist3->Scale(1./hist3->GetEntries());
  hist3->SetLineColor(kRed);
  hist3->SetMarkerStyle(21);
  hist3->SetMarkerColor(kRed);
  hs_3->Add(hist3);
  if(isLog)
    hs_3->SetMinimum(0.1/hist3->GetEntries());
  hs_3->Draw("ElpSAME");
  hs_3->GetXaxis()->SetTitle("# Displaced Tracks");
  hs_3->GetYaxis()->SetTitle("Normalized Collisions");

  TLegend *legend = new TLegend(0.75,0.75,0.90,0.90);
  legend->AddEntry(hist1,legendname_1.c_str(),"f");
  legend->AddEntry(hist2,(legendname_2+" not triggered").c_str(),"lep");
  legend->AddEntry(hist3,(legendname_2+" triggered").c_str(),"lep");
  legend->SetFillStyle(0);
  legend->Draw();


  pad2->cd();
  TH1F *hist_ratio2 = new TH1F("hist_ratio2","",nbins,min,max);
  hist_ratio2->SetMarkerStyle(21);
  hist_ratio2->SetMarkerColor(kBlue);
  hist_ratio2->SetLineColor(kBlue);
  hist_ratio2->GetXaxis()->SetTitle("# Displaced Tracks");
  hist_ratio2->GetYaxis()->SetTitle((legendname_2+"/"+legendname_1).c_str());
  hist_ratio2->GetXaxis()->SetLabelSize(.06);
  hist_ratio2->GetYaxis()->SetLabelSize(.06);
  hist_ratio2->GetXaxis()->SetTitleSize(.08);
  hist_ratio2->GetYaxis()->SetTitleSize(.08);
  hist_ratio2->GetYaxis()->SetTitleOffset(.55);
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
     hist_ratio2->SetBinContent(i,diff);
     hist_ratio2->SetBinError(i,diff_e);
   }
  hist_ratio2->Draw("ep");
  hist_ratio2->SetMinimum(0);
  hist_ratio2->SetMaximum(std::max(2.,hist_ratio2->GetMaximum()+hist_ratio2->GetBinError(hist_ratio2->GetMaximumBin())+0.1));

  TH1F *hist_ratio3 = new TH1F("hist_ratio3","",nbins,min,max);
  hist_ratio3->SetMarkerStyle(21);
  hist_ratio3->SetMarkerColor(kRed);
  hist_ratio3->SetLineColor(kRed);
  hist_ratio3->GetXaxis()->SetTitle("# Displaced Tracks");
  hist_ratio3->GetYaxis()->SetTitle((legendname_2+"/"+legendname_1).c_str());
  hist_ratio3->GetXaxis()->SetLabelSize(.06);
  hist_ratio3->GetYaxis()->SetLabelSize(.06);
  hist_ratio3->GetXaxis()->SetTitleSize(.08);
  hist_ratio3->GetYaxis()->SetTitleSize(.08);
  hist_ratio3->GetYaxis()->SetTitleOffset(.55);
   for (Int_t i=0;i<nbins;i++) {
     Double_t diff;
     Double_t diff_e;
     if (hist1->GetBinContent(i)){
       diff = hist3->GetBinContent(i)/hist1->GetBinContent(i);
       diff_e = diff*sqrt(pow(hist1->GetBinError(i)/hist1->GetBinContent(i),2)+pow(hist3->GetBinError(i)/hist3->GetBinContent(i),2));
     }
     else{
       diff = 0; diff_e = 0;
     }
     hist_ratio3->SetBinContent(i,diff);
     hist_ratio3->SetBinError(i,diff_e);
   }
  hist_ratio3->Draw("epSAME");
  hist_ratio2->SetMinimum(0);
  hist_ratio2->SetMaximum(std::max(2.,hist_ratio3->GetMaximum()+hist_ratio3->GetBinError(hist_ratio3->GetMaximumBin())+0.1));

  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kPink);
  line->Draw("SAME");
  pad2->Update();
  pad1->Update();
  c->Update();
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
