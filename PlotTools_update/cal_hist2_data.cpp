#include "helper.cpp"
void  cal_hist2_data(string dataname, string NVtx_cut,string fileout, string Obs1="sph", string Obs2="masspertrk", string massregion = "inclu", string nPU_cut= "inclu"){
  if(Obs1!="pu"&&Obs1!="sph"&&Obs1!="thrust"&&Obs1!="ndisplace"&&Obs1!="ntrk"&&Obs1!="pt"&&Obs1!="masspertrk"){
    cout<<"Invalid Observable"<<endl;
    return;
  }
  if(Obs2!="pu"&&Obs2!="sph"&&Obs2!="thrust"&&Obs2!="ndisplace"&&Obs2!="ntrk"&&Obs2!="pt"&&Obs2!="masspertrk"){
    cout<<"Invalid Observable"<<endl;
    return;
  }
  if(massregion!="inclu"&&massregion!="verylow"&&massregion!="verylowfull"&&massregion!="low"&&massregion!="lowfull"&&massregion!="medium"&&massregion!="high"&&massregion!="veryhigh"){
    cout<<"Invalide mass region"<<endl;
    return;
  }

  int nbins1=NBINS[Obs1];
  float min1=MIN[Obs1];
  float max1=MAX[Obs1];

  int nbins2=NBINS[Obs2];
  float min2=MIN[Obs2];
  float max2=MAX[Obs2];

  TChain *tdata = new TChain("ntuplizer/tree");
  tdata->AddFile(dataname.c_str());
//  tdata->AddFriend("tree",PUWeight_tree.c_str());
  tdata->AddFriend("tree",NVtx_cut.c_str());
//  if(MCPU!="") tdata->AddFriend("tree",MCPU.c_str());

  vector<int> *Obs1_int = new vector<int>[tdata->GetMaximum("PV_N")];
  vector<float> *Obs1_float = new vector<float>[tdata->GetMaximum("PV_N")];
  vector<int> *Obs2_int = new vector<int>[tdata->GetMaximum("PV_N")];
  vector<float> *Obs2_float = new vector<float>[tdata->GetMaximum("PV_N")];

  if(Obs1=="ntrk"||Obs1=="ndisplace") tdata->SetBranchAddress(BRANCHNAMEMC[Obs1].c_str(),&Obs1_int);
  else if(Obs1=="masspertrk"){
    tdata->SetBranchAddress(BRANCHNAMEMC["ntrk"].c_str(),&Obs1_int);
    tdata->SetBranchAddress("Instanton_Trk_TrkCut_mass",&Obs1_float);
  }
  else tdata->SetBranchAddress(BRANCHNAMEMC[Obs1].c_str(),&Obs1_float);
  if(Obs2=="ntrk"||Obs2=="ndisplace") tdata->SetBranchAddress(BRANCHNAMEMC[Obs2].c_str(),&Obs2_int);
  else if(Obs2=="masspertrk"){
    tdata->SetBranchAddress(BRANCHNAMEMC["ntrk"].c_str(),&Obs2_int);
    tdata->SetBranchAddress("Instanton_Trk_TrkCut_mass",&Obs2_float);
  }
  else tdata->SetBranchAddress(BRANCHNAMEMC[Obs2].c_str(),&Obs2_float);

//  float  PUWeight_N_PV_inclu_cut;
//  tdata->SetBranchAddress("PUWeight_N_PV_inclu_cut",&PUWeight_N_PV_inclu_cut);
  vector<bool> *is_inmasscut = new vector<bool>[tdata->GetMaximum("PV_N")];
  tdata->SetBranchAddress((MASSCUT[massregion]).c_str(),&is_inmasscut);
  int PV_N;
  tdata->SetBranchAddress("PV_N",&PV_N);

  TH2F *hist = new TH2F(("hist_"+Obs1+"_"+Obs2+"_"+massregion).c_str(),("Distribution of "+FULLNAME[Obs1]+" and "+FULLNAME[Obs2]+" in "+massregion+" mass").c_str(),nbins1,min1,max1,nbins2,min2,max2);
  hist->Sumw2(kTRUE);
  TH1F *hist1 = new TH1F("hist1",("Distribution of "+FULLNAME[Obs1]).c_str(),nbins1,min1,max1);
  hist1->Sumw2(kTRUE);
  TH1F *hist2 = new TH1F("hist2",("Distribution of "+FULLNAME[Obs2]).c_str(),nbins2,min2,max2);
  hist2->Sumw2(kTRUE);
  float Obs1_value,Obs2_value;
  float EX=0, EY=0, EXY=0, EX2=0, EY2=0;
  for(int i=0; i<tdata->GetEntries(); i++){
//  for(int i=0; i<1000; i++){
    tdata->GetEntry(i);
//    cout<<"Entry "<<i<<", PV_N "<<PV_N<<", cut len "<<is_inmasscut[0].size()<<endl;
    for(int iV=0; iV<PV_N; iV++){
      if(is_inmasscut[0][iV]){
        if(Obs1=="masspertrk") {
          if(Obs2=="ntrk") Obs1_value=Obs1_float[0][iV]/Obs2_int[0][iV];
          else Obs1_value=Obs1_float[0][iV]/Obs1_int[0][iV];
        }
        else if(Obs1=="ntrk"||Obs1=="ndisplace") Obs1_value=Obs1_int[0][iV];
        else Obs1_value=Obs1_float[0][iV];
        if(Obs2=="masspertrk") Obs2_value=Obs2_float[0][iV]/Obs2_int[0][iV];
        else if(Obs2=="ntrk"||Obs2=="ndisplace") Obs2_value=Obs2_int[0][iV];
        else Obs2_value=Obs2_float[0][iV];
//        cout<<"Obs1 "<<Obs1_value<<", Obs2 "<<Obs2_value<<", weight "<<PUWeight_N_PV_inclu_cut<<endl; 
        if(!isfinite(Obs1_value)||!isfinite(Obs2_value)) continue;
        hist->Fill(Obs1_value,Obs2_value);
        hist1->Fill(Obs1_value);
        hist2->Fill(Obs2_value);
        EX += Obs1_value;
        EY += Obs2_value;
        EXY += Obs1_value*Obs2_value;
        EX2 += Obs1_value*Obs1_value;
        EY2 += Obs2_value*Obs2_value;
      }
    }
  }
  EXY/=hist1->GetSumOfWeights();
  EX/=hist1->GetSumOfWeights();
  EY/=hist1->GetSumOfWeights();
  EX2/=hist1->GetSumOfWeights();
  EY2/=hist1->GetSumOfWeights();
  for(int ix=0; ix<=nbins1+1; ix++)
//    cout<<"hist1 bin "<<ix<<" "<<hist1->GetBinContent(ix)<<endl;
  for(int iy=0; iy<=nbins2+1; iy++)
//    cout<<"hist2 bin "<<iy<<" "<<hist2->GetBinContent(iy)<<endl;
//  cout<<"hist sum "<<hist->GetSumOfWeights()<<endl;
//  cout<<"hist1 sum "<<hist1->GetSumOfWeights()<<endl;
//  cout<<"hist2 sum "<<hist2->GetSumOfWeights()<<endl;
  hist->Scale(1./hist->GetSumOfWeights());
  hist1->Scale(1./hist1->GetSumOfWeights());
  hist2->Scale(1./hist2->GetSumOfWeights());
  double Pearson_corr = 0;
  Pearson_corr = (EXY-EX*EY)/sqrt((EX2-EX*EX)*(EY2-EY*EY));
  double corr;
  corr=hist->GetCorrelationFactor();
  double H_XY=0;
  double I_XY=0;
  double I_norm;
  for(int ix=1; ix<=nbins1; ix++){
    for(int iy=1; iy<=nbins2; iy++){
//      cout<<"ix = "<<ix<<", iy = "<<iy<<"hist "<<hist->GetBinContent(ix,iy)<<", hist1 "<<hist1->GetBinContent(ix)<<", hist2 "<<hist2->GetBinContent(iy)<<endl;
      if(hist->GetBinContent(ix,iy)>0){
        H_XY += -hist->GetBinContent(ix,iy)*log2(hist->GetBinContent(ix,iy));
        I_XY += hist->GetBinContent(ix,iy)*(log2(hist->GetBinContent(ix,iy))-log2(hist1->GetBinContent(ix))-log2(hist2->GetBinContent(iy)));
      }
    }
  }
  I_norm = I_XY/H_XY;
  
  TFile *f_out=TFile::Open(fileout.c_str(),"update");
  hist->Write();
  vector<double> result;
  result.push_back(corr);
  result.push_back(I_norm);
  result.push_back(Pearson_corr);
  f_out->WriteObject(&result,("corr_Inorm_"+Obs1+"_"+Obs2+"_"+massregion).c_str());
  f_out->Close();
  cout<<"Correlation "<<corr<<endl;
  cout<<"I_norm "<<I_norm<<endl;
  cout<<"Pearson_corr"<<Pearson_corr<<endl;
}
