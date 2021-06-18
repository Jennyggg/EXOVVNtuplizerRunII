#include <iostream>
#include <sstream>
#include <string>
#include <map>
using namespace std;
void splitstring(string origin, vector<string>& split ){
  istringstream iss(origin);
  string word;
  while (iss >> word){
    cout<<word<<"\n";
    split.push_back(word);
  }
  cout<<endl;
}
void rootfilelist(string origin, vector<string>& list){
  vector<string> split;
  splitstring(origin,split);
  for(int i=0; i<split.size(); i++){
    if(split[i].find(".root")==std::string::npos){
      FILE *proc = popen(("/bin/ls -1 "+split[i]+"/*.root").c_str(),"r");
      char buf[1024];
      while ( !feof(proc) && fgets(buf,sizeof(buf),proc) )
      {
//        printf("Line read: %s",buf);
          std::string s(buf);
          list.push_back(s.substr(0,s.size()-1));
          cout<<s<<endl;
      }
    }
    else{
      list.push_back(split[i]);
      cout<<split[i]<<endl;
    }
  }
}

void HistError(TChain * tree,string PUWeight_tree, string Branch_v, string Branch_error,string selection, TH1F * Hist_origin, int nbins, float min, float max){
  tree->AddFriend("tree",PUWeight_tree.c_str());
  TH1F *hist = new TH1F("hist","hist",nbins,min,max);
  hist->Sumw2(kTRUE);
  tree->Draw((Branch_v+">>hist").c_str(),("("+selection+")*"+Branch_error).c_str(),"goff");
  for (Int_t i=0;i<=nbins;i++) {
    Hist_origin->SetBinError(i,sqrt(Hist_origin->GetBinError(i)*Hist_origin->GetBinError(i)+hist->GetBinError(i)*hist->GetBinError(i)));
  }
}

map<string,int> NBINS{
{"sph",20},
{"thrust",20},
{"ndisplace",20},
{"ntrk",85},
{"ntrk_pt05",85},
{"ntrk_pt1",85},
{"pt",40},
{"pt_pt05",40},
{"pt_pt1",40},
{"masspertrk",35},
{"pu",65},
{"mass",50},
{"nheavyq",8},
{"ngendisplace",50},
{"ngenchargedisplace",50},
{"ngen",40},
{"ngencharged",35}
};
map<string,float> MIN{
{"sph",0},
{"thrust",0},
{"ndisplace",0},
{"ntrk",0},
{"ntrk_pt05",0},
{"ntrk_pt1",0},
{"pt",0},
{"pt_pt05",0},
{"pt_pt1",0},
{"masspertrk",0},
{"pu",0},
{"mass",0},
{"nheavyq",0},
{"ngendisplace",0},
{"ngenchargedisplace",0},
{"ngen",0},
{"ngencharged",0}
};
map<string,float> MAX{
{"sph",1},
{"thrust",0.5},
{"ndisplace",20},
{"ntrk",170},
{"ntrk_pt05",170},
{"ntrk_pt1",170},
{"pt",400},
{"pt_pt05",400},
{"pt_pt1",400},
{"masspertrk",35},
{"pu",65},
{"mass",1000},
{"nheavyq",8},
{"ngendisplace",50},
{"ngenchargedisplace",50},
{"ngen",400},
{"ngencharged",350}
};
map<string,string> BRANCHNAME{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"ntrk_pt05","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt05"},
{"ntrk_pt1","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt1"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"pt_pt05","Instanton_Trk_TrkCut_pt05_pt"},
{"pt_pt1","Instanton_Trk_TrkCut_pt1_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"mass","Instanton_Trk_TrkCut_mass"},
{"nheavyq","Instanton_N_genPromptHeavyQ"}
};
map<string,string> BRANCHNAMEMC{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
//{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
//{"ndisplace","Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity"},
{"ndisplace","Instanton_N_Trk_Displaced_fromPV2_highPurity"},
//{"ndisplace","Instanton_N_Trk_goodDisplaced_fromPV2_highPurity_sig3"},
//{"ndisplace","Instanton_N_Trk_fromPV2_highPurity_pt1"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"ntrk_pt05","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt05"},
{"ntrk_pt1","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt1"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"pt_pt05","Instanton_Trk_TrkCut_pt05_pt"},
{"pt_pt1","Instanton_Trk_TrkCut_pt1_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pu","nPuVtxTrue_"},
{"mass","Instanton_Trk_TrkCut_mass"},
{"nheavyq","Instanton_N_genPromptHeavyQ"},
{"ngendisplace","Instanton_N_gen_Displaced"},
{"ngenchargedisplace","Instanton_N_gen_ChargedDisplaced"},
{"ngen","Instanton_N_gen_FSParticle"},
{"ngencharged","Instanton_N_gen_ChargedFSParticle"}
};
map<string,string> BRANCHNAMEDATA{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"ntrk_pt05","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt05"},
{"ntrk_pt1","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt1"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"pt_pt05","Instanton_Trk_TrkCut_pt05_pt"},
{"pt_pt1","Instanton_Trk_TrkCut_pt1_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pu","Instan_Lumi_per_bunch_mean*69200"},
{"mass","Instanton_Trk_TrkCut_mass"}
};
map<string,string> PLOTNAME{
{"sph","Distribution of Spherocity"},
{"thrust","Distribution of Thrust"},
{"ndisplace","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity PVAssociationQuality>=4)"},
{"ntrk","Distribution of tracks (highPurity PVAssociationQuality>=4)"},
{"ntrk_pt05","Distribution of tracks (highPurity PVAssociationQuality>=4 pt>0.5 GeV)"},
{"ntrk_pt1","Distribution of tracks (highPurity PVAssociationQuality>=4 pt>1 GeV)"},
{"pt","Distribution of sum of track PT"},
{"pt_pt05","Distribution of sum of track PT (pt>0.5 GeV)"},
{"pt_pt1","Distribution of sum of track PT (pt>1 GeV)"},
{"masspertrk","Distribution of Mass/NTrack"},
{"pu","Distribution of Expected Pileup"},
{"mass","Distribution of C.O.M. mass of tracks"},
{"nheavyq","Distribution of #heavy quarks in the genevent"},
{"ngendisplace","Distribution of #displaced genparticles"},
{"ngenchargedisplace","Distribution of #charged displaced genparticles"},
{"ngen","Distribution of #genparticles"},
{"ngencharged","Distributio of #charged genparticles"}
};
map<string,string> AXISNAME{
{"sph","Spherocity"},
{"thrust","Thrust"},
{"ndisplace","NDisplace"},
{"ntrk","#Tracks"},
{"ntrk_pt05","#Tracks"},
{"ntrk_pt1","#Tracks"},
{"pt","sum of track PT (GeV)"},
{"pt_pt05","sum of track PT (GeV)"},
{"pt_pt1","sum of track PT (GeV)"},
{"masspertrk","Mass/NTrack (GeV)"},
{"pu","Pileup"},
{"mass","mass (GeV)"},
{"nheavyq","# Heavy quarks"},
{"ngendisplace","# displaced genparticle"},
{"ngenchargedisplace","# charged displaced genparticle"},
{"ngen","# genparticle"},
{"ngencharged","# charged genparticle"}
};
map<string,string> FULLNAME{
{"sph","Spherocity"},
{"thrust","Thrust"},
{"ndisplace","NDisplace"},
{"ntrk","#Tracks"},
{"ntrk_pt05","#Tracks"},
{"ntrk_pt1","#Tracks"},
{"pt","sum of track PT"},
{"pt_pt05","sum of track PT"},
{"pt_pt1","sum of track PT"},
{"masspertrk","Mass/NTrack"},
{"mass","Mass"},
{"nheavyq","# Heavy quarks"},
{"ngendisplace","# displaced genparticle"},
{"ngenchargedisplace","# charged displaced genparticle"},
{"ngen","# genparticle"},
{"ngencharged","# charged genparticle"}
};
map<string,string> MASSCUT{
{"inclu","PV_isinclu"},
{"verylow","PV_isverylowCR"},
{"verylowfull","PV_isverylow"},
{"low","PV_islowCR"},
{"lowfull","PV_islow"},
{"medium","PV_ismedium"},
{"high","PV_ishigh"},
{"veryhigh","PV_isveryhigh"}
};
map<string,string> LEGENDNAME{
{"inclu","           Inclusive mass region (full)"},
{"incluCR","           Inclusive mass region (CR)"},
{"verylow","           Very low mass region (CR)"},
{"verylowfull","           Very low mass region (full)"},
{"low","           Low mass region (CR)"},
{"lowfull","           Low mass region (full)"},
{"medium","           Medium mass region"},
{"high","           High mass region"},
{"veryhigh","         Very High mass region"}
};
map<string,string> SELE_MASS{
{"inclu","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0"},
{"incluCR","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&((Instanton_Trk_TrkCut_mass>=80)||(Instanton_Trk_TrkCut_mass<80&&Instanton_Trk_TrkCut_spherocity<0.3))"},
{"verylow","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&Instanton_Trk_TrkCut_spherocity<0.3"},
{"verylowfull","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40"},
{"low","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>40&&Instanton_Trk_TrkCut_mass<=80&&Instanton_Trk_TrkCut_spherocity<0.3"},
{"lowfull","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>40&&Instanton_Trk_TrkCut_mass<=80"},
{"medium","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>200&&Instanton_Trk_TrkCut_mass<=300"},
{"high","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0"},
{"veryhigh","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>300&&Instanton_Trk_TrkCut_mass<=500"}
};

