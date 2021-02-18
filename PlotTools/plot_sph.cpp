void plot_sph(string filename, string plotname, string legendname, bool isLog = false){
  int nbins=20;
  float min=0.;
  float max=1.;
  TChain *t = new TChain("ntuplizer/tree");
  t->AddFile(filename.c_str());
  TCanvas *c = new TCanvas("c","c",600,400);
  gStyle->SetOptStat(0);
  THStack *hs = new THStack("hs","Distribution of Spherocity");
  if(isLog){
    c->SetLogy();
    hs->SetMinimum(0.1);
  }
  TH1F *hist1 = new TH1F("hist1","Distribution of Spherocity",nbins,min,max);
//  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t->Draw("Instanton_Trk_spherocity>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB
  hist1=(TH1F*)gDirectory->Get("hist1");
  hist1->SetFillColor(kOrange);
  hist1->SetLineColor(kBlack);
  hs->Add(hist1);
  hs->Draw("hist");
  hs->Draw("ESAME");
  hs->GetXaxis()->SetTitle("Spherocity");
  hs->GetYaxis()->SetTitle("Events");
  TLegend *legend = new TLegend(0.75,0.80,0.90,0.90);
  legend->AddEntry(hist1,legendname.c_str(),"f");
  legend->SetFillStyle(0);
  legend->Draw();

  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
