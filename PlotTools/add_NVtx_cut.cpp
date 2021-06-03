void add_NVtx_cut(string filename, string NVtxName){
  TFile *f = new TFile(filename.c_str(),"READ");
  TTree *t = (TTree*)f->Get("ntuplizer/tree");
  TFile *f_v = new TFile(NVtxName.c_str(),"RECREATE");
  TTree *t_v = new TTree("tree","NVtx_pass_cut");
  int N_PV_inclu_cut, N_PV_verylow_cut, N_PV_low_cut, N_PV_medium_cut, N_PV_high_cut, N_PV_veryhigh_cut;
  t_v->Branch("N_PV_inclu_cut",&N_PV_inclu_cut,"N_PV_inclu_cut/I");
  t_v->Branch("N_PV_verylow_cut",&N_PV_verylow_cut,"N_PV_verylow_cut/I");
  t_v->Branch("N_PV_low_cut",&N_PV_low_cut,"N_PV_low_cut/I");
  t_v->Branch("N_PV_medium_cut",&N_PV_medium_cut,"N_PV_medium_cut/I");
  t_v->Branch("N_PV_high_cut",&N_PV_high_cut,"N_PV_high_cut/I");
  t_v->Branch("N_PV_veryhigh_cut",&N_PV_veryhigh_cut,"N_PV_veryhigh_cut/I");

  int PV_N,PV_N_good;
  vector<int> *Instanton_vtx_N_goodMuon = new vector<int>[t->GetMaximum("PV_N")];
  vector<bool> *PV_isgood = new vector<bool>[t->GetMaximum("PV_N")];
  vector<float> *Instanton_vtx_goodMuonIP = new vector<float>[t->GetMaximum("PV_N")];
  vector<float> *Instanton_Trk_TrkCut_mass = new vector<float>[t->GetMaximum("PV_N")];
  vector<float> *Instanton_Trk_TrkCut_spherocity = new vector<float>[t->GetMaximum("PV_N")];
  vector<float> *Instanton_vtx_genvertex_PVDistance = new vector<float>[t->GetMaximum("PV_N")];

  t->SetBranchAddress("Instanton_vtx_N_goodMuon",&Instanton_vtx_N_goodMuon);
  t->SetBranchAddress("PV_N",&PV_N);
  t->SetBranchAddress("PV_N_good",&PV_N_good);
  t->SetBranchAddress("PV_isgood",&PV_isgood);
  t->SetBranchAddress("Instanton_vtx_goodMuonIP",&Instanton_vtx_goodMuonIP);
  t->SetBranchAddress("Instanton_Trk_TrkCut_mass",&Instanton_Trk_TrkCut_mass);
  t->SetBranchAddress("Instanton_Trk_TrkCut_spherocity",&Instanton_Trk_TrkCut_spherocity);
//  t->SetBranchAddress("Instanton_vtx_genvertex_PVDistance",&Instanton_vtx_genvertex_PVDistance);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    N_PV_inclu_cut=N_PV_verylow_cut=N_PV_low_cut=N_PV_medium_cut=N_PV_high_cut=N_PV_veryhigh_cut=0;
    for(int iV=0; iV<PV_N; iV++){
//      if(!Instanton_vtx_N_goodMuon[0][iV]&&PV_isgood[0][iV]&&Instanton_vtx_goodMuonIP[0][iV]>2&&Instanton_Trk_TrkCut_mass[0][iV]>0&&Instanton_vtx_genvertex_PVDistance[0][iV]>2){
      if(!Instanton_vtx_N_goodMuon[0][iV]&&PV_isgood[0][iV]&&Instanton_vtx_goodMuonIP[0][iV]>2&&Instanton_Trk_TrkCut_mass[0][iV]>0){
        N_PV_inclu_cut++;
        if(Instanton_Trk_TrkCut_mass[0][iV]>20&&Instanton_Trk_TrkCut_mass[0][iV]<=40&&Instanton_Trk_TrkCut_spherocity[0][iV]<0.3) N_PV_verylow_cut++;
        else if(Instanton_Trk_TrkCut_mass[0][iV]>40&&Instanton_Trk_TrkCut_mass[0][iV]<=80&&Instanton_Trk_TrkCut_spherocity[0][iV]<0.3) N_PV_low_cut++;
        else if(Instanton_Trk_TrkCut_mass[0][iV]>80&&Instanton_Trk_TrkCut_mass[0][iV]<=200) N_PV_medium_cut++;
        else if(Instanton_Trk_TrkCut_mass[0][iV]>200&&Instanton_Trk_TrkCut_mass[0][iV]<=300) N_PV_high_cut++;
        else if(Instanton_Trk_TrkCut_mass[0][iV]>300&&Instanton_Trk_TrkCut_mass[0][iV]<=500) N_PV_veryhigh_cut++;
      }
    }
    t_v->Fill();
  }
  t_v->Write();
  f_v->Close();
}
