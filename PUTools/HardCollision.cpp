void HardCollision(string origin, string slim){
  TChain *t = new TChain("ntuplizer/tree");
  t->AddFile(origin.c_str());  
  t->SetBranchStatus("*", 1);
  for (auto disactiveBranchName : {"Instanton_N_Trk", "Instanton_N_Trk_Displaced", "Instanton_N_TrackJet","Instanton_Trk_mass","Instanton_Trk_spherocity","Instanton_Trk_broaden","Instanton_Trk_thrust"})
    t->SetBranchStatus(disactiveBranchName, 0);

  TFile *copy = new TFile(slim.c_str(),"recreate");
  TDirectory *dir = copy->mkdir("ntuplizer");
  dir->cd();
  
  TTree * t_copy = (TTree *) t->CloneTree();

  vector<vector<int>> Instanton_N_Trk;
  vector<vector<int>> Instanton_N_Trk_Displaced;
  vector<vector<int>> Instanton_N_TrackJet;
  vector<vector<float>> Instanton_Trk_mass;
  vector<vector<float>> Instanton_Trk_spherocity;
  vector<vector<float>> Instanton_Trk_broaden;
  vector<vector<float>> Instanton_Trk_thrust;

  for (auto activeBranchName : {"Instanton_N_Trk", "Instanton_N_Trk_Displaced", "Instanton_N_TrackJet","Instanton_Trk_mass","Instanton_Trk_spherocity","Instanton_Trk_broaden","Instanton_Trk_thrust"})
    t->SetBranchStatus(activeBranchName, 1);

  t->SetBranchAddress("Instanton_N_Trk",(vector<int> *)Instanton_N_Trk);
  t->SetBranchAddress("Instanton_N_Trk_Displaced",(vector<int> *)Instanton_N_Trk_Displaced);
  t->SetBranchAddress("Instanton_N_TrackJet",(vector<int> *)Instanton_N_TrackJet);
  t->SetBranchAddress("Instanton_Trk_mass",(vector<float> *)Instanton_Trk_mass);
  t->SetBranchAddress("Instanton_Trk_spherocity",(vector<float> *)Instanton_Trk_spherocity);
  t->SetBranchAddress("Instanton_Trk_broaden",(vector<float> *)Instanton_Trk_broaden);
  t->SetBranchAddress("Instanton_Trk_thrust",(vector<float> *)Instanton_Trk_thrust);

  int Instanton_N_Trk_s;
  int Instanton_N_Trk_Displaced_s;
  int Instanton_N_TrackJet_s;
  float Instanton_Trk_mass_s;
  float Instanton_Trk_spherocity_s;
  float Instanton_Trk_broaden_s;
  float Instanton_Trk_thrust_s;

  TBranch * Instanton_N_Trk_b = t_copy->Branch("Instanton_N_Trk",&Instanton_N_Trk_s,"Instanton_N_Trk/I");
  TBranch * Instanton_N_Trk_Displaced_b = t_copy->Branch("Instanton_N_Trk_Displaced",&Instanton_N_Trk_Displaced_s,"Instanton_N_Trk_Displaced/I");
  TBranch * Instanton_N_TrackJet_b = t_copy->Branch("Instanton_N_TrackJet",&Instanton_N_TrackJet_s,"Instanton_N_TrackJet/I");
  TBranch * Instanton_Trk_mass_b = t_copy->Branch("Instanton_Trk_mass",&Instanton_Trk_mass_s,"Instanton_Trk_mass/F");
  TBranch * Instanton_Trk_spherocity_b = t_copy->Branch("Instanton_Trk_spherocity",&Instanton_Trk_spherocity_s,"Instanton_Trk_spherocity/F");
  TBranch * Instanton_Trk_broaden_b = t_copy->Branch("Instanton_Trk_broaden",&Instanton_Trk_broaden_s,"Instanton_Trk_broaden/F");
  TBranch * Instanton_Trk_thrust_b = t_copy->Branch("Instanton_Trk_thrust",&Instanton_Trk_thrust_s,"Instanton_Trk_thrust/F");

  for(long long int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    Instanton_N_Trk_s = Instanton_N_Trk[0][0];
    Instanton_N_Trk_Displaced_s = Instanton_N_Trk_Displaced[0][0];
    Instanton_N_TrackJet_s = Instanton_N_TrackJet[0][0];
    Instanton_Trk_mass_s = Instanton_Trk_mass[0][0];
    Instanton_Trk_spherocity_s = Instanton_Trk_spherocity[0][0];
    Instanton_Trk_broaden_s = Instanton_Trk_broaden[0][0];
    Instanton_Trk_thrust_s = Instanton_Trk_thrust[0][0];
    Instanton_N_Trk_b->Fill();
    Instanton_N_Trk_Displaced_b->Fill();
    Instanton_N_TrackJet_b->Fill();
    Instanton_Trk_mass_b->Fill();
    Instanton_Trk_spherocity_b->Fill();
    Instanton_Trk_broaden_b->Fill();
    Instanton_Trk_thrust_b->Fill();
  }
  copy->Write();
}

