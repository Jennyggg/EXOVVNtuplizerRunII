void add_NVtx_cut(string filename, string NVtxName, bool isvtxcut, bool isPU){
  TFile *f = new TFile(filename.c_str(),"READ");
  TTree *t = (TTree*)f->Get("ntuplizer/tree");
  TFile *f_v = new TFile(NVtxName.c_str(),"RECREATE");
  TTree *t_v = new TTree("tree","NVtx_pass_cut");
  int N_PV_inclu_cut, N_PV_verylow_cut, N_PV_low_cut, N_PV_medium_cut, N_PV_high_cut, N_PV_veryhigh_cut;
  vector<bool> PV_isinclu, PV_isverylow, PV_islow, PV_isverylowCR, PV_islowCR, PV_ismedium, PV_ishigh, PV_isveryhigh;

  t_v->Branch("N_PV_inclu_cut",&N_PV_inclu_cut,"N_PV_inclu_cut/I");
  t_v->Branch("N_PV_verylow_cut",&N_PV_verylow_cut,"N_PV_verylow_cut/I");
  t_v->Branch("N_PV_low_cut",&N_PV_low_cut,"N_PV_low_cut/I");
  t_v->Branch("N_PV_medium_cut",&N_PV_medium_cut,"N_PV_medium_cut/I");
  t_v->Branch("N_PV_high_cut",&N_PV_high_cut,"N_PV_high_cut/I");
  t_v->Branch("N_PV_veryhigh_cut",&N_PV_veryhigh_cut,"N_PV_veryhigh_cut/I");

  t_v->Branch("PV_isinclu",&PV_isinclu);
  t_v->Branch("PV_isverylow",&PV_isverylow);
  t_v->Branch("PV_islow",&PV_islow);
  t_v->Branch("PV_isverylowCR",&PV_isverylowCR);
  t_v->Branch("PV_islowCR",&PV_islowCR);
  t_v->Branch("PV_ismedium",&PV_ismedium);
  t_v->Branch("PV_ishigh",&PV_ishigh);
  t_v->Branch("PV_isveryhigh",&PV_isveryhigh);

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
  if(isvtxcut)
    t->SetBranchAddress("Instanton_vtx_genvertex_PVDistance",&Instanton_vtx_genvertex_PVDistance);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    N_PV_inclu_cut=N_PV_verylow_cut=N_PV_low_cut=N_PV_medium_cut=N_PV_high_cut=N_PV_veryhigh_cut=0;
    float mindis_genvtx=999.;
    int iV_mindis_genvtx=0;
    if(isvtxcut && !isPU){  // to choose the reconstructed hardest collision (the primary vertex) in minbias MC, select the vertex closest to the genvertex
      for(int iV=0; iV<PV_N; iV++){
        if(Instanton_vtx_genvertex_PVDistance[0][iV]<mindis_genvtx){
          mindis_genvtx=Instanton_vtx_genvertex_PVDistance[0][iV];
          iV_mindis_genvtx=iV;
        }
      }
    }
    for(int iV=0; iV<PV_N; iV++){

      PV_isinclu.push_back(false);
      PV_isverylow.push_back(false);
      PV_islow.push_back(false);
      PV_isverylowCR.push_back(false);
      PV_islowCR.push_back(false);
      PV_ismedium.push_back(false);
      PV_ishigh.push_back(false);
      PV_isveryhigh.push_back(false);
//      if(!Instanton_vtx_N_goodMuon[0][iV]&&PV_isgood[0][iV]&&Instanton_vtx_goodMuonIP[0][iV]>2&&Instanton_Trk_TrkCut_mass[0][iV]>0&&Instanton_vtx_genvertex_PVDistance[0][iV]>2){
      if(!Instanton_vtx_N_goodMuon[0][iV]&&PV_isgood[0][iV]&&Instanton_vtx_goodMuonIP[0][iV]>2&&Instanton_Trk_TrkCut_mass[0][iV]>0){
        if(isvtxcut && isPU && Instanton_vtx_genvertex_PVDistance[0][iV]<=2) continue; 
        if(isvtxcut && !isPU && iV!=iV_mindis_genvtx) continue;
        N_PV_inclu_cut++;
        PV_isinclu[iV]=true;
        if(Instanton_Trk_TrkCut_mass[0][iV]>20&&Instanton_Trk_TrkCut_mass[0][iV]<=40){
          PV_isverylow[iV]=true;
          if(Instanton_Trk_TrkCut_spherocity[0][iV]<0.3) {
            N_PV_verylow_cut++;
            PV_isverylowCR[iV]=true;
          }
        }
        else if(Instanton_Trk_TrkCut_mass[0][iV]>40&&Instanton_Trk_TrkCut_mass[0][iV]<=80){

          PV_islow[iV]=true;
          if(Instanton_Trk_TrkCut_spherocity[0][iV]<0.3) {
            N_PV_low_cut++;
            PV_islowCR[iV]=true;
          }
        }
        else if(Instanton_Trk_TrkCut_mass[0][iV]>80&&Instanton_Trk_TrkCut_mass[0][iV]<=200){ 
          N_PV_medium_cut++;
          PV_ismedium[iV]=true;
        }
        else if(Instanton_Trk_TrkCut_mass[0][iV]>200&&Instanton_Trk_TrkCut_mass[0][iV]<=300){ 
          N_PV_high_cut++;
          PV_ishigh[iV]=true;
        }
        else if(Instanton_Trk_TrkCut_mass[0][iV]>300&&Instanton_Trk_TrkCut_mass[0][iV]<=500){ 
          N_PV_veryhigh_cut++;
          PV_isveryhigh[iV]=true;
        }
      }
    }
    t_v->Fill();
    PV_isinclu.clear();
    PV_isverylow.clear();
    PV_isverylowCR.clear();
    PV_islow.clear();
    PV_islowCR.clear();
    PV_ismedium.clear();
    PV_ishigh.clear();
    PV_isveryhigh.clear();
  }
  t_v->Write();
  f_v->Close();
}
