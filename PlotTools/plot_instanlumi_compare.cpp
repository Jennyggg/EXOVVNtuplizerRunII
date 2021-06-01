#include "helper.cpp"
void plot_instanlumi_compare(string filenames_1,string filenames_2, string plotname, string legendname_1,string legendname_2, bool isLog = false){
  int nbins=60;
  float min=0.00015;
  float max=0.00075;
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
  THStack *hs_1 = new THStack("hs_1","Distribution of Instantaneous Luminosity per bunch");
  THStack *hs_2 = new THStack("hs_2","Distribution of Instantaneous Luminosity per bunch");
//  THStack *hs_1 = new THStack("hs_1","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity)");
//  THStack *hs_2 = new THStack("hs_2","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity)");
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
  TH1F *hist = new TH1F("hist","Distribution of Instantaneous Luminosity per bunch",nbins,min,max);
  TH1F *hist1 = new TH1F("hist1","Distribution of Instantaneous Luminosity per bunch",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of Instantaneous Luminosity per bunch",nbins,min,max);

  string selection = "PV_N_good>0";
  for(int i=0; i<filename_1.size(); i++){
    t->AddFile(filename_1[i].c_str());
    if(t->GetEntries(selection.c_str())>0){
      cout<<"processing file "<<filename_1[i]<<endl;
      t->Draw("Instan_Lumi_per_bunch_mean>>hist",selection.c_str(),"goff");

      hist=(TH1F*)gDirectory->Get("hist");
      hist1->Add(hist);
    }
    t->Reset();
  }
  for(int i=0; i<filename_2.size(); i++){
    t->AddFile(filename_2[i].c_str());
    if(t->GetEntries(selection.c_str())>0){
      cout<<"processing file "<<filename_2[i]<<endl;
      t->Draw("Instan_Lumi_per_bunch_mean>>hist",selection.c_str(),"goff");
      hist=(TH1F*)gDirectory->Get("hist");
      hist2->Add(hist);
    }
    t->Reset();
  }


  hist1->Scale(1./hist1->GetSumOfWeights());
  hist1->SetFillColor(kOrange);
//  hist1->SetLineColor(kBlack);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("Instantaneous Luminosity per bunch (#mub^{-1})");
  hs_1->GetYaxis()->SetTitle("Normalized Events");

  hist2->Scale(1./hist2->GetSumOfWeights());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("Instantaneous Luminosity per bunch (#mub^{-1})");
  hs_2->GetYaxis()->SetTitle("Normalized Events");
  hs_1->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));
  hs_2->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist2,"UUNORM");
  TLegend *legend = new TLegend(0.70,0.8,0.90,0.90);
  legend->AddEntry(hist1,legendname_1.c_str(),"f");
  legend->AddEntry(hist2,legendname_2.c_str(),"lep");
//  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
//  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();


  pad2->cd();
  TH1F *hist_ratio = new TH1F("hist_ratio","",nbins,min,max);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->SetMarkerColor(kBlack);
  hist_ratio->GetXaxis()->SetTitle("Instantaneous Luminosity per bunch (#mub^{-1})");
  hist_ratio->GetYaxis()->SetTitle((legendname_2+"/"+legendname_1).c_str());
  hist_ratio->GetXaxis()->SetLabelSize(.06);
  hist_ratio->GetYaxis()->SetLabelSize(.06);
  hist_ratio->GetXaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleOffset(.55);
   for (Int_t i=0;i<=nbins;i++) {
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
  pad2->Update();
  pad1->Update();
  c->Update();
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
