void plot_IPparticle(string filename="../Ntuplizer/flatTuple.root"){
  TChain *tMC = new TChain("ntuplizer/tree");
  tMC->AddFile(filename.c_str());
  TCanvas *c = new TCanvas("c","c",800,600);
  c->SetLogy();
//  THStack *hs1 = new THStack("hs1","Distribution of transverse IP");
//  THStack *hs2 = new THStack("hs2","Distribution of transverse IP");
//  THStack *hs3 = new THStack("hs3","Distribution of transverse IP");
//  TH1F *hist1 = new TH1F("hist1","Distribution of transverse IP",30,0,0.3);
//  TH1F *hist2 = new TH1F("hist2","Distribution of transverse IP",30,0,0.3);
//  TH1F *hist3 = new TH1F("hist3","Distribution of transverse IP",30,0,0.3);

//  THStack *hs1 = new THStack("hs1","Distribution of longitudinal IP");
//  THStack *hs2 = new THStack("hs2","Distribution of longitudinal IP");
//  THStack *hs3 = new THStack("hs3","Distribution of longitudinal IP");
//  TH1F *hist1 = new TH1F("hist1","Distribution of longitudinal IP",40,0,40);
//  TH1F *hist2 = new TH1F("hist2","Distribution of longitudinal IP",40,0,40);
//  TH1F *hist3 = new TH1F("hist3","Distribution of longitudinal IP",40,0,40);

  THStack *hs1 = new THStack("hs1","Distribution of 3D IP");
  THStack *hs2 = new THStack("hs2","Distribution of 3D IP");
  THStack *hs3 = new THStack("hs3","Distribution of 3D IP");
  TH1F *hist1 = new TH1F("hist1","Distribution of 3D IP",40,0,1.);
  TH1F *hist2 = new TH1F("hist2","Distribution of 3D IP",40,0,1.);
  TH1F *hist3 = new TH1F("hist3","Distribution of 3D IP",40,0,1.);


//  tMC->Draw("genParticle_dxy>>hist1","genParticle_isfromheavyq==1","goff");
//  tMC->Draw("genParticle_dxy>>hist2","genParticle_isfromheavyq==0","goff");
//  tMC->Draw("genParticle_dxy>>hist3","","goff");
//  tMC->Draw("genParticle_dz>>hist1","genParticle_isfromheavyq==1","goff");
//  tMC->Draw("genParticle_dz>>hist2","genParticle_isfromheavyq==0","goff");
//  tMC->Draw("genParticle_dz>>hist3","","goff");
  tMC->Draw("genParticle_IP3D>>hist1","genParticle_isfromheavyq==1","goff");
  tMC->Draw("genParticle_IP3D>>hist2","genParticle_isfromheavyq==0","goff");
  tMC->Draw("genParticle_IP3D>>hist3","","goff");

  hist1=(TH1F*)gDirectory->Get("hist1");
  hist2=(TH1F*)gDirectory->Get("hist2");
  hist3=(TH1F*)gDirectory->Get("hist3");
  hist1->Scale(1./hist1->GetSumOfWeights());
  hist2->Scale(1./hist2->GetSumOfWeights());
  hist3->Scale(1./hist3->GetSumOfWeights());
  hs1->SetMinimum(0.1/hist1->GetEntries());
  hs2->SetMinimum(0.1/hist2->GetEntries());
  hs3->SetMinimum(0.1/hist3->GetEntries());

  hist1->SetLineColor(kRed);
  hist1->SetMarkerStyle(21);
  hist1->SetMarkerColor(kRed);
  hs1->Add(hist1);
  hist2->SetLineColor(kBlue);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);
  hs2->Add(hist2);
  hist3->SetFillColor(kOrange);
  hs3->Add(hist3);

  hs3->Draw("hist");
  hs3->Draw("ESAME");
//  hs3->GetXaxis()->SetTitle("transversse IP (cm)");
//  hs3->GetXaxis()->SetTitle("longitudinal IP (cm)");
  hs3->GetXaxis()->SetTitle("3D IP (cm)");
  hs3->GetYaxis()->SetTitle("Normalized particles");

  hs1->Draw("ElpSAME");
  hs2->Draw("ElpSAME");

  hs3->SetMaximum(hist2->GetMaximum());

  TLegend *legend = new TLegend(0.65,0.65,0.90,0.90);
  legend->AddEntry(hist3,"all particles","f");
  legend->AddEntry(hist1,"from heavy quarks","lep");
  legend->AddEntry(hist2,"from light partons","lep");
  legend->Draw();

  c->Update();
//  c->SaveAs("plots/MinBias_transIPgenparticle.png");
//  c->SaveAs("plots/MinBias_transIPgenparticle.pdf");
//  c->SaveAs("plots/MinBias_longIPgenparticle.png");
//  c->SaveAs("plots/MinBias_longIPgenparticle.pdf");
  c->SaveAs("plots/MinBias_3dIPgenparticle.png");
  c->SaveAs("plots/MinBias_3dIPgenparticle.pdf");
}
