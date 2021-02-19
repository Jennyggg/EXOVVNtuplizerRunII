#include "helper.cpp"
void plot_thrust_compare(string filenames_1,string filenames_2, string plotname, string legendname_1,string legendname_2, bool isLog = false){
  int nbins=20;
  float min=0.;
  float max=0.5;
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
  TCanvas *c = new TCanvas("c","c",600,400);
  gStyle->SetOptStat(0);
  THStack *hs_1 = new THStack("hs_1","Distribution of Thrust");
  THStack *hs_2 = new THStack("hs_2","Distribution of Thrust");
  if(isLog){
    c->SetLogy();
//    hs_1->SetMinimum(0.1);
//    hs_2->SetMinimum(0.1);
  }
  TH1F *hist1 = new TH1F("hist1","Distribution of Thrust",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of Thrust",nbins,min,max);
  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_1->Draw("Instanton_Trk_thrust>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_2->Draw("Instanton_Trk_thrust>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

  hist1=(TH1F*)gDirectory->Get("hist1");
  hist1->Scale(1./hist1->GetEntries());
  hist1->SetFillColor(kOrange);
//  hist1->SetLineColor(kBlack);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("Thrust");
  hs_1->GetYaxis()->SetTitle("Normalized Events");

  hist2=(TH1F*)gDirectory->Get("hist2");
  hist2->Scale(1./hist2->GetEntries());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("Thrust");
  hs_2->GetYaxis()->SetTitle("Normalized Events");
  TLegend *legend = new TLegend(0.75,0.80,0.90,0.90);
  legend->AddEntry(hist1,legendname_1.c_str(),"f");
  legend->AddEntry(hist2,legendname_2.c_str(),"lep");
  legend->SetFillStyle(0);
  legend->Draw();

  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
