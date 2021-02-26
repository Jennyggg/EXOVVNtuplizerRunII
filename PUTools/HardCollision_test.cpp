void HardCollision_test(string origin, string slim){
  TChain *t = new TChain("ntuplizer/tree");
  t->AddFile(origin.c_str());  

  TFile *copy = new TFile(slim.c_str(),"recreate");
  TDirectory *dir = copy->mkdir("ntuplizer");
  dir->cd();
  
  TTree * t_copy = new TTree("tree","tree");

  vector<int> *Instanton_N_Trk = new vector<int>;
  vector<int> *Instanton_N_Trk_Displaced = new vector<int>;
  vector<int> *Instanton_N_TrackJet = new vector<int>;
  vector<float> *Instanton_Trk_mass = new vector<float>;
  vector<float> *Instanton_Trk_spherocity = new vector<float>;
  vector<float> *Instanton_Trk_broaden = new vector<float>;
  vector<float> *Instanton_Trk_thrust = new vector<float>;


  t->SetBranchAddress("Instanton_N_Trk",&Instanton_N_Trk);
  t->SetBranchAddress("Instanton_N_Trk_Displaced",&Instanton_N_Trk_Displaced);
  t->SetBranchAddress("Instanton_N_TrackJet",&Instanton_N_TrackJet);
  t->SetBranchAddress("Instanton_Trk_mass",&Instanton_Trk_mass);
  t->SetBranchAddress("Instanton_Trk_spherocity",&Instanton_Trk_spherocity);
  t->SetBranchAddress("Instanton_Trk_broaden",&Instanton_Trk_broaden);
  t->SetBranchAddress("Instanton_Trk_thrust",&Instanton_Trk_thrust);

  int Instanton_N_Trk_s;
  int Instanton_N_Trk_Displaced_s;
  int Instanton_N_TrackJet_s;
  float Instanton_Trk_mass_s;
  float Instanton_Trk_spherocity_s;
  float Instanton_Trk_broaden_s;
  float Instanton_Trk_thrust_s;

  TBranch * Instanton_N_Trk_b = t_copy->Branch("Instanton_N_Trk_hard",&Instanton_N_Trk_s,"Instanton_N_Trk_hard/I");
  TBranch * Instanton_N_Trk_Displaced_b = t_copy->Branch("Instanton_N_Trk_Displaced_hard",&Instanton_N_Trk_Displaced_s,"Instanton_N_Trk_Displaced_hard/I");
  TBranch * Instanton_N_TrackJet_b = t_copy->Branch("Instanton_N_TrackJet_hard",&Instanton_N_TrackJet_s,"Instanton_N_TrackJet_hard/I");
  TBranch * Instanton_Trk_mass_b = t_copy->Branch("Instanton_Trk_mass_hard",&Instanton_Trk_mass_s,"Instanton_Trk_mass_hard/F");
  TBranch * Instanton_Trk_spherocity_b = t_copy->Branch("Instanton_Trk_spherocity_hard",&Instanton_Trk_spherocity_s,"Instanton_Trk_spherocity_hard/F");
  TBranch * Instanton_Trk_broaden_b = t_copy->Branch("Instanton_Trk_broaden_hard",&Instanton_Trk_broaden_s,"Instanton_Trk_broaden_hard/F");
  TBranch * Instanton_Trk_thrust_b = t_copy->Branch("Instanton_Trk_thrust_hard",&Instanton_Trk_thrust_s,"Instanton_Trk_thrust_hard/F");

  for(long long int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    Instanton_N_Trk_s = (*Instanton_N_Trk)[0];
    Instanton_N_Trk_Displaced_s = (*Instanton_N_Trk_Displaced)[0];
    Instanton_N_TrackJet_s = (*Instanton_N_TrackJet)[0];
    Instanton_Trk_mass_s = (*Instanton_Trk_mass)[0];
    Instanton_Trk_spherocity_s = (*Instanton_Trk_spherocity)[0];
    Instanton_Trk_broaden_s = (*Instanton_Trk_broaden)[0];
    Instanton_Trk_thrust_s = (*Instanton_Trk_thrust)[0];
    t_copy->Fill();
//    Instanton_N_Trk_b->Fill();
//    Instanton_N_Trk_Displaced_b->Fill();
//    Instanton_N_TrackJet_b->Fill();
//    Instanton_Trk_mass_b->Fill();
//    Instanton_Trk_spherocity_b->Fill();
//   Instanton_Trk_broaden_b->Fill();
//    Instanton_Trk_thrust_b->Fill();
//    cout<<"Filled"<<endl;
  }
  t_copy->Write("", TObject::kOverwrite);
  dir->Write();
  copy->Write();
}

