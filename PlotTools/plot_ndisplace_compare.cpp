#include "helper.cpp"
void plot_ndisplace_compare(string filenames_1,string filenames_2, string plotname, string legendname_1,string legendname_2, bool isLog = false){
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
//  t_1->AddFriend("ntuplizer/tree","/pnfs/psi.ch/cms/trivcat/store/user/jinw/ZeroBias2018/flatTuple_small_add_Iso.root");
//  t_2->AddFriend("ntuplizer/tree","/pnfs/psi.ch/cms/trivcat/store/user/jinw/ZeroBias2018/flatTuple_small_add_Iso.root");
  TCanvas *c = new TCanvas("c","c",600,600);
  gStyle->SetOptStat(0);
  THStack *hs_1 = new THStack("hs_1","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity PVAssociationQuality>=4)");
  THStack *hs_2 = new THStack("hs_2","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity PVAssociationQuality>=4)");
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
  TH1F *hist1 = new TH1F("hist1","Distribution of #displaced tracks",nbins,min,max);
  TH1F *hist2 = new TH1F("hist2","Distribution of #displaced tracks",nbins,min,max);

//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","(PV_N_good)==15&&Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","(PV_N_good)==20&&Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&(PV_N_good)==30&&Instanton_vtx_goodMuonIP>2&&EVENT_lumiBlock>=550&&EVENT_lumiBlock<600&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","Instanton_vtx_goodMuonIP>2&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","EVENT_lumiBlock>=550&&EVENT_lumiBlock<600&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","Instanton_vtx_goodMuonIP>2&&EVENT_lumiBlock>=1490&&EVENT_lumiBlock<1507&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_highPurity>>hist1","Instanton_vtx_goodMuonIP>2&&EVENT_lumiBlock>=25&&EVENT_lumiBlock<65&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","EVENT_lumiBlock>=25&&EVENT_lumiBlock<65&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","PV_N==25&&EVENT_lumiBlock>=25&&EVENT_lumiBlock<65&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist1","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("Instanton_N_Trk_goodDisplaced_highPurity>>hist1","EVENT_run==316569&&EVENT_lumiBlock<600&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");//to compare the behavior of different luminosities
//  t_1->Draw("(Instanton_N_Trk_goodDisplaced_PVAssociationQuality4+Instanton_N_Trk_goodDisplaced_PVAssociationQuality5+Instanton_N_Trk_goodDisplaced_PVAssociationQuality6+Instanton_N_Trk_goodDisplaced_PVAssociationQuality7)>>hist1","PV_Iso>0.4&&EVENT_run==316569&&EVENT_lumiBlock<600&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("(Instanton_N_Trk_goodDisplaced)>>hist1","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_1->Draw("(Instanton_N_Trk_Displaced_PVAssociationQuality4+Instanton_N_Trk_Displaced_PVAssociationQuality5+Instanton_N_Trk_Displaced_PVAssociationQuality6+Instanton_N_Trk_Displaced_PVAssociationQuality7)>>hist1","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); 
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","PV_N>20&&PV_N<25&&!Instanton_vtx_N_goodMuon&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_1->Draw("Instanton_N_Trk_Displaced>>hist1","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB


//    t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood");

//    t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","(PV_N_good)==15&&Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood");
//    t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","(PV_N_good)==20&&Instanton_vtx_goodMuonIP>2&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//    t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","(PV_N_good)==20&&Instanton_vtx_goodMuonIP>2&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&(PV_N_good)==30&&Instanton_vtx_goodMuonIP>2&&EVENT_lumiBlock>=550&&EVENT_lumiBlock<600&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","EVENT_lumiBlock>=550&&EVENT_lumiBlock<600&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","Instanton_vtx_goodMuonIP>2&&EVENT_lumiBlock>=1490&&EVENT_lumiBlock<1507&&EVENT_run==316569&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","EVENT_lumiBlock>=25&&EVENT_lumiBlock<65&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","PV_N==25&&EVENT_lumiBlock>=25&&EVENT_lumiBlock<65&&EVENT_run==315974&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity>>hist2","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("Instanton_N_Trk_goodDisplaced_highPurity>>hist2","EVENT_run==316569&&EVENT_lumiBlock>1300&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("(Instanton_N_Trk_goodDisplaced_PVAssociationQuality4+Instanton_N_Trk_goodDisplaced_PVAssociationQuality5+Instanton_N_Trk_goodDisplaced_PVAssociationQuality6+Instanton_N_Trk_goodDisplaced_PVAssociationQuality7)>>hist2","EVENT_run==316569&&EVENT_lumiBlock>1300&&!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("(Instanton_N_Trk_goodDisplaced)>>hist2","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40");
//  t_2->Draw("(Instanton_N_Trk_Displaced_PVAssociationQuality4+Instanton_N_Trk_Displaced_PVAssociationQuality5+Instanton_N_Trk_Displaced_PVAssociationQuality6+Instanton_N_Trk_Displaced_PVAssociationQuality7)>>hist2","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); 
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","PV_N>20&&PV_N<25&&!Instanton_vtx_N_goodMuon&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total
//  t_2->Draw("Instanton_N_Trk_goodDisplaced>>hist2","!Instanton_vtx_isBPHtrigger&&Instanton_Trk_mass>20&&Instanton_Trk_mass<=40"); //total exclude BPH trigger collisions
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0"); //standard
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //event shape
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced>6"); //tight
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>15&&(Instanton_Trk_mass/Instanton_N_Trk)>2.0&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3"); //CRA
//  t_2->Draw("Instanton_N_Trk_Displaced>>hist2","Instanton_Trk_mass>20&&Instanton_Trk_mass<=40&&Instanton_N_Trk>20&&(Instanton_Trk_mass/Instanton_N_Trk)<1.5&&Instanton_N_TrackJet==0&&Instanton_Trk_thrust>0.3&&Instanton_N_Trk_Displaced<4"); //CRB

  hist1=(TH1F*)gDirectory->Get("hist1");
  hist1->Scale(1./hist1->GetSumOfWeights());
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
  hist2->Scale(1./hist2->GetSumOfWeights());
  hist2->SetLineColor(kRed);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs_2->Add(hist2);
  if(isLog)
    hs_2->SetMinimum(0.1/hist2->GetEntries());
  hs_2->Draw("ElpSAME");
  hs_2->GetXaxis()->SetTitle("# Displaced Tracks");
  hs_2->GetYaxis()->SetTitle("Normalized Collisions");

  Double_t chi2_ndof = 0;
  Double_t p_value = 0;
  chi2_ndof = hist1->Chi2Test(hist2,"UUNORMCHI2/NDF");
  p_value = hist1->Chi2Test(hist2,"UUNORM");
  TLegend *legend = new TLegend(0.70,0.65,0.90,0.90);
  legend->AddEntry(hist1,legendname_1.c_str(),"f");
  legend->AddEntry(hist2,legendname_2.c_str(),"lep");
  legend->AddEntry((TObject*)0, ("chi2/ndof = "+ std::to_string(chi2_ndof)).c_str(), "");
  legend->AddEntry((TObject*)0, ("p-value = "+ std::to_string(p_value)).c_str(), "");
  legend->SetFillStyle(0);
  legend->Draw();


  pad2->cd();
  TH1F *hist_ratio = new TH1F("hist_ratio","",nbins,min,max);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->SetMarkerColor(kBlack);
  hist_ratio->GetXaxis()->SetTitle("# Displaced Tracks");
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
