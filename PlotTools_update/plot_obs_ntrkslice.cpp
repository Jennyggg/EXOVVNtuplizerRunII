#include "helper.cpp"
void plot_obs_ntrkslice(string resultname, string Obs, string massregion, string plotname,bool isLog = false){
  TFile * f = TFile::Open(resultname.c_str(),"READ");
  TH2F * hist;
  int nbins=NBINS[Obs];
  float min=MIN[Obs];
  float max=MAX[Obs];
  TH1F * hist0 = new TH1F(("hist_"+Obs+"_ntrk20-40_"+massregion).c_str(),("Distribution of "+Obs+" with ntrk 20-40").c_str(),nbins,min,max);
  TH1F * hist1 = new TH1F(("hist_"+Obs+"_ntrk40-60_"+massregion).c_str(),("Distribution of "+Obs+" with ntrk 40-60").c_str(),nbins,min,max);
  TH1F * hist2 = new TH1F(("hist_"+Obs+"_ntrk60-80_"+massregion).c_str(),("Distribution of "+Obs+" with ntrk 60-80").c_str(),nbins,min,max);
  TH1F * hist3 = new TH1F(("hist_"+Obs+"_ntrk80-100_"+massregion).c_str(),("Distribution of "+Obs+" with ntrk 80-100").c_str(),nbins,min,max);
  TH1F * hist4 = new TH1F(("hist_"+Obs+"_ntrk100-120_"+massregion).c_str(),("Distribution of "+Obs+" with ntrk 100-120").c_str(),nbins,min,max);

//  if(Obs=="sph"||Obs=="thrust"||Obs=="ndisplace"||Obs=="pt_pt05"){
  if(Obs=="sph"||Obs=="thrust"||Obs=="ndisplace"||Obs=="pt"){
    f->GetObject(("hist_ntrk_"+Obs+"_"+massregion).c_str(),hist);
//    f->GetObject(("hist_ntrk_pt05_"+Obs+"_"+massregion).c_str(),hist);
//    f->GetObject(("hist_ntrk_pt1_"+Obs+"_"+massregion).c_str(),hist);
    for(int ibin=0; ibin <= nbins+1; ibin++){
      float hist0_iibin_value=0;
      float hist0_iibin_error=0;
      for (int iibin=5; iibin <= 8; iibin++){
        hist0_iibin_value += hist->GetBinContent(iibin,ibin);
        hist0_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist0_iibin_error = sqrt(hist0_iibin_error);
      hist0->SetBinContent(ibin,hist0_iibin_value);
      hist0->SetBinError(ibin,hist0_iibin_error);

      float hist1_iibin_value=0;
      float hist1_iibin_error=0;
      for (int iibin=9; iibin <= 12; iibin++){
        hist1_iibin_value += hist->GetBinContent(iibin,ibin);
        hist1_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist1_iibin_error = sqrt(hist1_iibin_error);
      hist1->SetBinContent(ibin,hist1_iibin_value);
      hist1->SetBinError(ibin,hist1_iibin_error);

      float hist2_iibin_value=0;
      float hist2_iibin_error=0;
      for (int iibin=13; iibin <= 16; iibin++){
        hist2_iibin_value += hist->GetBinContent(iibin,ibin);
        hist2_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist2_iibin_error = sqrt(hist2_iibin_error);
      hist2->SetBinContent(ibin,hist2_iibin_value);
      hist2->SetBinError(ibin,hist2_iibin_error);

      float hist3_iibin_value=0;
      float hist3_iibin_error=0;
      for (int iibin=17; iibin <= 20; iibin++){
        hist3_iibin_value += hist->GetBinContent(iibin,ibin);
        hist3_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist3_iibin_error = sqrt(hist3_iibin_error);
      hist3->SetBinContent(ibin,hist3_iibin_value);
      hist3->SetBinError(ibin,hist3_iibin_error);

      float hist4_iibin_value=0;
      float hist4_iibin_error=0;
      for (int iibin=21; iibin <= 24; iibin++){
        hist4_iibin_value += hist->GetBinContent(iibin,ibin);
        hist4_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist4_iibin_error = sqrt(hist4_iibin_error);
      hist4->SetBinContent(ibin,hist4_iibin_value);
      hist4->SetBinError(ibin,hist4_iibin_error);
    }
    hist0->Scale(1./hist0->GetSumOfWeights());
    hist1->Scale(1./hist1->GetSumOfWeights());
    hist2->Scale(1./hist2->GetSumOfWeights());
    hist3->Scale(1./hist3->GetSumOfWeights());
    hist4->Scale(1./hist4->GetSumOfWeights());
  }
  else{
    f->GetObject(("hist_"+Obs+"_ntrk_"+massregion).c_str(),hist);
//    f->GetObject(("hist_"+Obs+"_ntrk_pt05_"+massregion).c_str(),hist);
//    f->GetObject(("hist_"+Obs+"_ntrk_pt1_"+massregion).c_str(),hist);
    for(int ibin=0; ibin <= nbins+1; ibin++){
      float hist0_iibin_value=0;
      float hist0_iibin_error=0;
      for (int iibin=5; iibin <= 8; iibin++){
        hist0_iibin_value += hist->GetBinContent(ibin,iibin);
        hist0_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist0_iibin_error = sqrt(hist0_iibin_error);
      hist0->SetBinContent(ibin,hist0_iibin_value);
      hist0->SetBinError(ibin,hist0_iibin_error);

      float hist1_iibin_value=0;
      float hist1_iibin_error=0;
      for (int iibin=9; iibin <= 12; iibin++){
        hist1_iibin_value += hist->GetBinContent(ibin,iibin);
        hist1_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist1_iibin_error = sqrt(hist1_iibin_error);
      hist1->SetBinContent(ibin,hist1_iibin_value);
      hist1->SetBinError(ibin,hist1_iibin_error);

      float hist2_iibin_value=0;
      float hist2_iibin_error=0;
      for (int iibin=13; iibin <= 16; iibin++){
        hist2_iibin_value += hist->GetBinContent(ibin,iibin);
        hist2_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist2_iibin_error = sqrt(hist2_iibin_error);
      hist2->SetBinContent(ibin,hist2_iibin_value);
      hist2->SetBinError(ibin,hist2_iibin_error);

      float hist3_iibin_value=0;
      float hist3_iibin_error=0;
      for (int iibin=17; iibin <= 20; iibin++){
        hist3_iibin_value += hist->GetBinContent(ibin,iibin);
        hist3_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist3_iibin_error = sqrt(hist3_iibin_error);
      hist3->SetBinContent(ibin,hist3_iibin_value);
      hist3->SetBinError(ibin,hist3_iibin_error);

      float hist4_iibin_value=0;
      float hist4_iibin_error=0;
      for (int iibin=21; iibin <= 24; iibin++){
        hist4_iibin_value += hist->GetBinContent(ibin,iibin);
        hist4_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist4_iibin_error = sqrt(hist4_iibin_error);
      hist4->SetBinContent(ibin,hist4_iibin_value);
      hist4->SetBinError(ibin,hist4_iibin_error);
    }
    hist0->Scale(1./hist0->GetSumOfWeights());
    hist1->Scale(1./hist1->GetSumOfWeights());
    hist2->Scale(1./hist2->GetSumOfWeights());
    hist3->Scale(1./hist3->GetSumOfWeights());
    hist4->Scale(1./hist4->GetSumOfWeights());
  }
  hist0->SetFillColor(kOrange);
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  THStack *hs_0 = new THStack("hs_0",PLOTNAME[Obs].c_str());
  THStack *hs_1 = new THStack("hs_1",PLOTNAME[Obs].c_str());
  THStack *hs_2 = new THStack("hs_2",PLOTNAME[Obs].c_str());
  THStack *hs_3 = new THStack("hs_3",PLOTNAME[Obs].c_str());
  THStack *hs_4 = new THStack("hs_4",PLOTNAME[Obs].c_str());
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  if(isLog){
    pad1->SetLogy();
    hs_0->SetMinimum(0.1/hist->GetEntries());
    hs_1->SetMinimum(0.1/hist->GetEntries());
    hs_2->SetMinimum(0.1/hist->GetEntries());
    hs_3->SetMinimum(0.1/hist->GetEntries());
    hs_4->SetMinimum(0.1/hist->GetEntries());
  }
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  hs_0->Add(hist0);
  hs_0->Draw("hist");
  hs_0->Draw("ESAME");
  hs_0->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_0->GetYaxis()->SetTitle("Normalized Collisions");
  hist1->SetLineColor(kBlue);
  hist1->SetMarkerStyle(21);
  hist1->SetMarkerColor(kBlue);
  hs_1->Add(hist1);
  hs_1->Draw("ElpSAME");
  hs_1->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_1->GetYaxis()->SetTitle("Normalized Collisions");

  hist2->SetLineColor(kGreen);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kGreen);
  hs_2->Add(hist2);
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_2->GetYaxis()->SetTitle("Normalized Collisions");

  hist3->SetLineColor(kRed);
  hist3->SetMarkerStyle(21);
  hist3->SetMarkerColor(kRed);
  hs_3->Add(hist3);
  hs_3->Draw("ElpSAME");
  hs_3->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_3->GetYaxis()->SetTitle("Normalized Collisions");

  hist4->SetLineColor(kBlack);
  hist4->SetMarkerStyle(21);
  hist4->SetMarkerColor(kBlack);
  hs_4->Add(hist4);
  hs_4->Draw("ElpSAME");
  hs_4->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_4->GetYaxis()->SetTitle("Normalized Collisions");

  double y_max=std::max(hist0->GetMaximum(),hist1->GetMaximum());
  y_max= std::max(y_max, hist2->GetMaximum());
  y_max= std::max(y_max, hist3->GetMaximum());
  y_max= std::max(y_max, hist4->GetMaximum());
  hs_0->SetMaximum(y_max);
  hs_1->SetMaximum(y_max);
  hs_2->SetMaximum(y_max);
  hs_3->SetMaximum(y_max);
  hs_4->SetMaximum(y_max);


  Double_t chi2_ndof01 = 0;
  Double_t p_value01 = 0;
  chi2_ndof01 = hist1->Chi2Test(hist0,"UUNORMCHI2/NDF");
  p_value01 = hist1->Chi2Test(hist0,"UUNORM");
  Double_t chi2_ndof12 = 0;
  Double_t p_value12 = 0;
  chi2_ndof12 = hist2->Chi2Test(hist1,"UUNORMCHI2/NDF");
  p_value12 = hist2->Chi2Test(hist1,"UUNORM");
  Double_t chi2_ndof23 = 0;
  Double_t p_value23 = 0;
  chi2_ndof23 = hist3->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value23 = hist3->Chi2Test(hist2,"UUNORM");
  Double_t chi2_ndof34 = 0;
  Double_t p_value34 = 0;
  chi2_ndof34 = hist4->Chi2Test(hist3,"UUNORMCHI2/NDF");
  p_value34 = hist4->Chi2Test(hist3,"UUNORM");

  TLegend *legend = new TLegend(0.65,0.45,0.90,0.90);
  legend->AddEntry(hist0,"ntrack 20-40","f");
  legend->AddEntry(hist1,"ntrack 40-60","lep");
  legend->AddEntry(hist2,"ntrack 60-80","lep");
  legend->AddEntry(hist3,"ntrack 80-100","lep");
  legend->AddEntry(hist4,"ntrack 100-120","lep");
  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof01)+", " + std::to_string(chi2_ndof12)).c_str(), "");
  legend->AddEntry((TObject*)0, (std::to_string(chi2_ndof23)+", "+std::to_string(chi2_ndof34)).c_str(), "");
  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value01) + ", "+ std::to_string(p_value12)).c_str(), "");
  legend->AddEntry((TObject*)0, (std::to_string(p_value23)+", "+std::to_string(p_value34)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();
  pad2->cd();

  TH1F *hist_diff1 = new TH1F("hist_diff1","",nbins,min,max);
  hist_diff1->SetMarkerStyle(21);
  hist_diff1->SetMarkerColor(kBlue);
  hist_diff1->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist_diff1->GetYaxis()->SetTitle("ratio");
  hist_diff1->GetXaxis()->SetLabelSize(.06);
  hist_diff1->GetYaxis()->SetLabelSize(.06);
  hist_diff1->GetXaxis()->SetTitleSize(.08);
  hist_diff1->GetYaxis()->SetTitleSize(.08);
  hist_diff1->GetYaxis()->SetTitleOffset(.55);

  TH1F *hist_diff2 = new TH1F("hist_diff2","",nbins,min,max);
  hist_diff2->SetMarkerStyle(21);
  hist_diff2->SetMarkerColor(kGreen);
  hist_diff2->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist_diff2->GetYaxis()->SetTitle("ratio");
  hist_diff2->GetXaxis()->SetLabelSize(.06);
  hist_diff2->GetYaxis()->SetLabelSize(.06);
  hist_diff2->GetXaxis()->SetTitleSize(.08);
  hist_diff2->GetYaxis()->SetTitleSize(.08);
  hist_diff2->GetYaxis()->SetTitleOffset(.55);

  TH1F *hist_diff3 = new TH1F("hist_diff3","",nbins,min,max);
  hist_diff3->SetMarkerStyle(21);
  hist_diff3->SetMarkerColor(kRed);
  hist_diff3->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist_diff3->GetYaxis()->SetTitle("ratio");
  hist_diff3->GetXaxis()->SetLabelSize(.06);
  hist_diff3->GetYaxis()->SetLabelSize(.06);
  hist_diff3->GetXaxis()->SetTitleSize(.08);
  hist_diff3->GetYaxis()->SetTitleSize(.08);
  hist_diff3->GetYaxis()->SetTitleOffset(.55);

  TH1F *hist_diff4 = new TH1F("hist_diff4","",nbins,min,max);
  hist_diff4->SetMarkerStyle(21);
  hist_diff4->SetMarkerColor(kBlack);
  hist_diff4->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist_diff4->GetYaxis()->SetTitle("ratio");
  hist_diff4->GetXaxis()->SetLabelSize(.06);
  hist_diff4->GetYaxis()->SetLabelSize(.06);
  hist_diff4->GetXaxis()->SetTitleSize(.08);
  hist_diff4->GetYaxis()->SetTitleSize(.08);
  hist_diff4->GetYaxis()->SetTitleOffset(.55);

  for (Int_t i=0;i<nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    if (hist0->GetBinContent(i)&&hist1->GetBinContent(i)){
      diff = hist1->GetBinContent(i)/hist0->GetBinContent(i);
      diff_e = diff*sqrt(pow(hist0->GetBinError(i)/hist0->GetBinContent(i),2)+pow(hist1->GetBinError(i)/hist1->GetBinContent(i),2));
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist_diff1->SetBinContent(i,diff);
    hist_diff1->SetBinError(i,diff_e); 
  }

  for (Int_t i=0;i<nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    if (hist0->GetBinContent(i)&&hist2->GetBinContent(i)){
      diff = hist2->GetBinContent(i)/hist0->GetBinContent(i);
      diff_e = diff*sqrt(pow(hist0->GetBinError(i)/hist0->GetBinContent(i),2)+pow(hist2->GetBinError(i)/hist2->GetBinContent(i),2));
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist_diff2->SetBinContent(i,diff);
    hist_diff2->SetBinError(i,diff_e);
  }

  for (Int_t i=0;i<nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    if (hist0->GetBinContent(i)&&hist3->GetBinContent(i)){
      diff = hist3->GetBinContent(i)/hist0->GetBinContent(i);
      diff_e = diff*sqrt(pow(hist0->GetBinError(i)/hist0->GetBinContent(i),2)+pow(hist3->GetBinError(i)/hist3->GetBinContent(i),2));
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist_diff3->SetBinContent(i,diff);
    hist_diff3->SetBinError(i,diff_e);
  }

  for (Int_t i=0;i<nbins;i++) {
    Double_t diff;
    Double_t diff_e;
    if (hist0->GetBinContent(i)&&hist4->GetBinContent(i)){
      diff = hist4->GetBinContent(i)/hist0->GetBinContent(i);
      diff_e = diff*sqrt(pow(hist0->GetBinError(i)/hist0->GetBinContent(i),2)+pow(hist4->GetBinError(i)/hist4->GetBinContent(i),2));
    }
    else{
      diff = 0; diff_e = 0;
    }
    hist_diff4->SetBinContent(i,diff);
    hist_diff4->SetBinError(i,diff_e);
  }

  hist_diff1->Draw("p");
  hist_diff2->Draw("pSAME");
  hist_diff3->Draw("pSAME");
  hist_diff4->Draw("pSAME");
  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
  pad2->Update();
  pad1->Update();
  c->Update();
  plotname += "_"+massregion+"_"+Obs+"_ntrkslice";
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
  delete c;

}
