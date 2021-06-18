#include "helper.cpp"
void plot_obs_ndisplaceslice(string resultname, string Obs, string massregion, string plotname,bool isLog = false){
  TFile * f = TFile::Open(resultname.c_str(),"READ");
  TH2F * hist;
  int nbins=NBINS[Obs];
  float min=MIN[Obs];
  float max=MAX[Obs];
  TH1F * hist0 = new TH1F(("hist_"+Obs+"_ndisplace0_"+massregion).c_str(),("Distribution of "+Obs+" with ndisplace = 0").c_str(),nbins,min,max);
  TH1F * hist1 = new TH1F(("hist_"+Obs+"_ndisplacegt0_"+massregion).c_str(),("Distribution of "+Obs+" with ndisplace > 0").c_str(),nbins,min,max);
  if(Obs=="sph"||Obs=="thrust"){
    f->GetObject(("hist_ndisplace_"+Obs+"_"+massregion).c_str(),hist);
    for(int ibin=0; ibin <= nbins+1; ibin++){
      hist0->SetBinContent(ibin,hist->GetBinContent(1,ibin));
      hist0->SetBinError(ibin,hist->GetBinError(1,ibin));
      float hist1_iibin_value=0;
      float hist1_iibin_error=0;
      for (int iibin=2; iibin <= NBINS["ndisplace"]+1; iibin++){
        hist1_iibin_value += hist->GetBinContent(iibin,ibin);
        hist1_iibin_error += hist->GetBinError(iibin,ibin)*hist->GetBinError(iibin,ibin);
      }
      hist1_iibin_error = sqrt(hist1_iibin_error);
      hist1->SetBinContent(ibin,hist1_iibin_value);
      hist1->SetBinError(ibin,hist1_iibin_error);
    }
    hist0->Scale(1./hist0->GetSumOfWeights());
    hist1->Scale(1./hist1->GetSumOfWeights());
  }
  else{
    f->GetObject(("hist_"+Obs+"_ndisplace_"+massregion).c_str(),hist);
    for(int ibin=0; ibin <= nbins+1; ibin++){
      hist0->SetBinContent(ibin,hist->GetBinContent(ibin,1));
      hist0->SetBinError(ibin,hist->GetBinError(ibin,1));
      float hist1_iibin_value=0;
      float hist1_iibin_error=0;
      for (int iibin=2; iibin <= NBINS["ndisplace"]+1; iibin++){
        hist1_iibin_value += hist->GetBinContent(ibin,iibin);
        hist1_iibin_error += hist->GetBinError(ibin,iibin)*hist->GetBinError(ibin,iibin);
      }
      hist1_iibin_error = sqrt(hist1_iibin_error);
      hist1->SetBinContent(ibin,hist1_iibin_value);
      hist1->SetBinError(ibin,hist1_iibin_error);
    }
    hist0->Scale(1./hist0->GetSumOfWeights());
    hist1->Scale(1./hist1->GetSumOfWeights());

  }
  hist0->SetFillColor(kOrange);
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  THStack *hs_0 = new THStack("hs_0",PLOTNAME[Obs].c_str());
  THStack *hs_1 = new THStack("hs_1",PLOTNAME[Obs].c_str());
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  if(isLog){
    pad1->SetLogy();
    hs_0->SetMinimum(0.1/hist->GetEntries());
  }
  hs_0->Add(hist0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  hs_0->Draw("hist");
  hs_0->Draw("ESAME");
  hs_0->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_0->GetYaxis()->SetTitle("Normalized Collisions");
  hist1->SetLineColor(kRed);
  hist1->SetMarkerStyle(21);
  hist1->SetMarkerColor(kBlue);
  hs_1->Add(hist1);
  if(isLog)
    hs_1->SetMinimum(0.1/hist->GetEntries());
  hs_1->Draw("ElpSAME");
  hs_1->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hs_1->GetYaxis()->SetTitle("Normalized Collisions");
  hs_0->SetMaximum(std::max(hist0->GetMaximum(),hist1->GetMaximum()));
  hs_1->SetMaximum(std::max(hist0->GetMaximum(),hist1->GetMaximum()));

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist0,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist0,"UUNORM");

  TLegend *legend = new TLegend(0.65,0.65,0.90,0.90);
  legend->AddEntry(hist0,"ndisplace=0","f");
  legend->AddEntry(hist1,"ndisplace>0","lep");
  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();
  pad2->cd();
  TH1F *hist_diff = new TH1F("hist_diff","",nbins,min,max);
  hist_diff->SetMarkerStyle(21);
  hist_diff->SetMarkerColor(kBlack);
  hist_diff->GetXaxis()->SetTitle(AXISNAME[Obs].c_str());
  hist_diff->GetYaxis()->SetTitle("ratio");
  hist_diff->GetXaxis()->SetLabelSize(.06);
  hist_diff->GetYaxis()->SetLabelSize(.06);
  hist_diff->GetXaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleSize(.08);
  hist_diff->GetYaxis()->SetTitleOffset(.55);
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
    hist_diff->SetBinContent(i,diff);
    hist_diff->SetBinError(i,diff_e); 
  }
  hist_diff->Draw("p");
  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
  pad2->Update();
  pad1->Update();
  c->Update();
  plotname += "_"+massregion+"_"+Obs+"_ndisplaceslice";
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
  delete c;

}
