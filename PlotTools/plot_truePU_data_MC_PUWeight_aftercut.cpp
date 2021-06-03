#include "helper.cpp"
void plot_truePU_data_MC_PUWeight_aftercut(string MCname,string dataname, string PUWeight_tree, string MC_NVtx, string data_NVtx, string plotname, string MClegend,string datalegend,string selectionMC="", string selectiondata="",bool PUreweight = true, bool isLog = false){
  int nbins=65;
  float min=0;
  float max=65;
  TChain *tMC = new TChain("ntuplizer/tree");
  TChain *tdata = new TChain("ntuplizer/tree");
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  THStack *hs_1 = new THStack("hs_1","Distribution of pileup expectation"); //MC
  THStack *hs_2 = new THStack("hs_2","Distribution of average pileup in a lumiblock");            //data
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
  TH1F *hist = new TH1F("hist","Distribution of true pileup",nbins,min,max);
  TH1F *hist1 = new TH1F("hist1","Distribution of true pileup",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of true pileup",nbins,min,max);
  hist->Sumw2(kTRUE);
  tMC->AddFile(MCname.c_str());
  tMC->AddFriend("tree",PUWeight_tree.c_str());
  tMC->AddFriend("tree",MC_NVtx.c_str());
  tdata->AddFriend("tree",data_NVtx.c_str());
    if(tMC->GetEntries(selectionMC.c_str())>0){
      if(PUreweight){
 //       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*PV_N_good").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_inclu_cut").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_verylow_cut").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_low_cut").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_medium_cut").c_str(),"goff");
//       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_high_cut").c_str(),"goff");
       tMC->Draw("nPuVtxTrue>>hist",("("+selectionMC+")*PUWeight*N_PV_veryhigh_cut").c_str(),"goff");
        hist=(TH1F*)gDirectory->Get("hist");
        hist1->Add(hist);
        HistError(tMC,PUWeight_tree,"nPuVtxTrue","PUWeightE",selectionMC,hist1,nbins,min,max);
      }
      else{
        tMC->Draw("nPuVtxTrue>>hist",selectionMC.c_str(),"goff");
        hist=(TH1F*)gDirectory->Get("hist");
        hist1->Add(hist);
      }
    }

  tdata->AddFile(dataname.c_str());    
    if(tdata->GetEntries(selectiondata.c_str())>0){
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",selectiondata.c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*PV_N_good").c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_inclu_cut").c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_verylow_cut").c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_low_cut").c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_medium_cut").c_str(),"goff");
//      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_high_cut").c_str(),"goff");
      tdata->Draw("Instan_Lumi_per_bunch_mean*69200>>hist",("("+selectiondata+")*N_PV_veryhigh_cut").c_str(),"goff");
      hist=(TH1F*)gDirectory->Get("hist");
      hist2->Add(hist);
    }


//  hist1->Scale(1./hist1->GetSumOfWeights());
//  hist1->Scale(1./3.80443e+07);
  hist1->Scale(1./6.98185e+07);
  cout<<"hist 1 Sum of Weights"<<hist1->GetSumOfWeights()<<endl;
//  hist1->Scale(1./hist1->GetSumOfWeights());
  hist1->SetFillColor(kOrange);
//  hist1->SetLineColor(kBlack);
  if(isLog)
    hs_1->SetMinimum(0.1/hist1->GetEntries());
  hs_1->Add(hist1);
  hs_1->Draw("hist");
  hs_1->Draw("ESAME");
  hs_1->GetXaxis()->SetTitle("True Pileup");
  hs_1->GetYaxis()->SetTitle("Normalized Events");

//  hist2->Scale(1./hist2->GetSumOfWeights());
//  hist2->Scale(1./6.90581e+07);
  hist2->Scale(1./6.90581e+07);
  cout<<"hist 2 Sum of Weights"<<hist2->GetSumOfWeights()<<endl;
//  hist2->Scale(1./hist2->GetSumOfWeights());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("True Pileup");
  hs_2->GetYaxis()->SetTitle("Normalized Events");
  hs_1->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));
  hs_2->SetMaximum(std::max(hist1->GetMaximum(),hist2->GetMaximum()));

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist2,"UUNORM");
  TLegend *legend = new TLegend(0.70,0.8,0.90,0.90);
  legend->AddEntry(hist1,MClegend.c_str(),"f");
  legend->AddEntry(hist2,datalegend.c_str(),"lep");
//  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
//  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();


  pad2->cd();
  TH1F *hist_ratio = new TH1F("hist_ratio","",nbins,min,max);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->SetMarkerColor(kBlack);
  hist_ratio->GetXaxis()->SetTitle("True Pileup");
  hist_ratio->GetYaxis()->SetTitle((datalegend+"/"+MClegend).c_str());
  hist_ratio->GetXaxis()->SetLabelSize(.06);
  hist_ratio->GetYaxis()->SetLabelSize(.06);
  hist_ratio->GetXaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleSize(.08);
  hist_ratio->GetYaxis()->SetTitleOffset(.55);
   for (Int_t i=0;i<=nbins;i++) {
     Double_t diff;
     Double_t diff_e;
     if (hist1->GetBinContent(i)&&hist2->GetBinContent(i)){
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
//  hist_ratio->SetMaximum(std::max(2.,hist_ratio->GetMaximum()+hist_ratio->GetBinError(hist_ratio->GetMaximumBin())+0.1));
  hist_ratio->SetMaximum(2);
  TLine *line = new TLine(min,1,max,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
  pad2->Update();
  pad1->Update();
  c->Update();
  c->SaveAs((plotname+".png").c_str());
  c->SaveAs((plotname+".pdf").c_str());
}
