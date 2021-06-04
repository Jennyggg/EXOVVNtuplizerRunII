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
{"pt",40},
{"masspertrk",35},
{"pu",65}
};
map<string,float> MIN{
{"sph",0},
{"thrust",0},
{"ndisplace",0},
{"ntrk",0},
{"pt",0},
{"masspertrk",0},
{"pu",0}
};
map<string,float> MAX{
{"sph",1},
{"thrust",0.5},
{"ndisplace",20},
{"ntrk",170},
{"pt",400},
{"masspertrk",35},
{"pu",65}
};
map<string,string> BRANCHNAME{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"}
};
map<string,string> BRANCHNAMEMC{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pu","nPuVtxTrue_"}
};
map<string,string> BRANCHNAMEDATA{
{"sph","Instanton_Trk_TrkCut_spherocity"},
{"thrust","Instanton_Trk_TrkCut_thrust"},
{"ndisplace","Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity"},
{"ntrk","Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pt","Instanton_Trk_TrkCut_pt"},
{"masspertrk","Instanton_Trk_TrkCut_mass/Instanton_N_Trk_PVAssociationQualityLeq4_highPurity"},
{"pu","Instan_Lumi_per_bunch_mean*69200"}
};
map<string,string> PLOTNAME{
{"sph","Distribution of Spherocity"},
{"thrust","Distribution of Thrust"},
{"ndisplace","Distribution of #displaced tracks (IP>0.02cm IP/sigma(IP)>5 highPurity PVAssociationQuality>=4)"},
{"ntrk","Distribution of tracks (highPurity PVAssociationQuality>=4)"},
{"pt","Distribution of sum of track PT"},
{"masspertrk","Distribution of Mass/NTrack"}
};
map<string,string> AXISNAME{
{"sph","Spherocity"},
{"thrust","Thrust"},
{"ndisplace","NDisplace"},
{"ntrk","#Tracks"},
{"pt","sum of track PT (GeV)"},
{"masspertrk","Mass/NTrack (GeV)"}
};
map<string,string> LEGENDNAME{
{"inclu","           Inclusive mass region"},
{"verylow","           Very low mass region"},
{"low","           Low mass region"},
{"medium","           Medium mass region"},
{"high","           High mass region"},
{"veryhigh","         Very High mass region"}
};
map<string,string> SELE_MASS{
{"inclu","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0"},
{"verylow","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>20&&Instanton_Trk_TrkCut_mass<=40&&Instanton_Trk_TrkCut_spherocity<0.3"},
{"low","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>40&&Instanton_Trk_TrkCut_mass<=80&&Instanton_Trk_TrkCut_spherocity<0.3"},
{"medium","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>200&&Instanton_Trk_TrkCut_mass<=300"},
{"high","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0"},
{"veryhigh","!Instanton_vtx_N_goodMuon&&PV_isgood&&Instanton_vtx_goodMuonIP>2&&Instanton_Trk_TrkCut_mass>0&&Instanton_Trk_TrkCut_mass>300&&Instanton_Trk_TrkCut_mass<=500"}
};

