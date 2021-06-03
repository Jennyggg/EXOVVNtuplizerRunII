void correct_nPuVtxTrue(string filename, string NVtxName){
  TFile *f = new TFile(filename.c_str(),"READ");
  TTree *t = (TTree*)f->Get("ntuplizer/tree");
  TFile *f_v = new TFile(NVtxName.c_str(),"RECREATE");
  TTree *t_v = new TTree("tree","nPuVtxTrue");
  float nPuVtxTrue_;
  t_v->Branch("nPuVtxTrue_",&nPuVtxTrue_,"nPuVtxTrue_/F");
  vector<float> *nPuVtxTrue = new vector<float>[10];

  t->SetBranchAddress("nPuVtxTrue",&nPuVtxTrue);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    nPuVtxTrue_=nPuVtxTrue[0][0];
    t_v->Fill();
  } 
  t_v->Write();
  f_v->Close();
}
